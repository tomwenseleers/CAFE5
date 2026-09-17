#include "innovation.h"
#include "clade.h"
#include "optimizer.h"
#include "optimizer_scorer.h"
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <memory>
#include <numeric>
#include <random>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>
#ifdef _OPENMP
#include <omp.h>
#endif

namespace innovation {
using Vec = std::vector<double>;
const double INF = std::numeric_limits<double>::infinity();

std::vector<double> transition(int maximum, double lambda, double nu, double time) {
    if (maximum < 0 || !std::isfinite(lambda) || !std::isfinite(nu) ||
        !std::isfinite(time) || lambda < 0 || nu < 0 || time < 0)
        throw std::runtime_error("Invalid BDI transition argument");
    if (!std::isfinite(lambda*time) || !std::isfinite(nu*time))
        throw std::runtime_error("Rate-time product exceeds numerical range");
    const int s = maximum + 1;
    Vec a(size_t(s)*s, 0);
    if (time == 0) { for (int i=0;i<s;++i) a[size_t(i)*s+i]=1; return a; }
    if (lambda == 0) {
        Vec p(s); double lp=-nu*time; p[0]=std::exp(lp);
        for (int j=1;j<s;++j) {lp+=(nu>0?std::log(nu*time):-INF)-std::log(double(j));p[j]=std::exp(lp);}
        for (int i=0;i<s;++i) for (int j=i;j<s;++j) a[size_t(i)*s+j]=p[j-i];
        return a;
    }
    double x=lambda*time, q=x/(1+x), survival=1/(1+x);
    // Immigration PGF [1+lambda*t*(1-z)]^(-nu/lambda).
    // Avoid nu/lambda in recurrence, retaining the Poisson limit numerically.
    double lp=-nu*time*(x==0?1:std::log1p(x)/x); a[0]=std::exp(lp);
    for (int j=1;j<s;++j) {
        double numerator=nu*time+(j-1)*x;
        lp+=(numerator>0?std::log(numerator):-INF)-std::log(double(j))-std::log1p(x);
        a[j]=std::exp(lp);
    }
    // Multiply by the one-ancestor BD PGF for each successive starting count.
    // Geometric convolution accumulator makes this O(s^2), all terms positive.
    for (int i=1;i<s;++i) {
        double sum=0;
        for (int j=0;j<s;++j) {
            if (j) sum=a[size_t(i-1)*s+j-1]+q*sum;
            a[size_t(i)*s+j]=q*a[size_t(i-1)*s+j]+survival*survival*sum;
        }
    }
    return a;
}

struct Options {
    std::string tree,input,prefix="innovation",root_file,matrix_file;
    double lambda=-1,nu=-1,root_mean=-1,matrix_time=-1;
    int maximum=80,iterations=400,starts=3,simulate=0,bootstrap=0,threads=1;
    unsigned seed=20260917;
    bool unconditioned=false,check=true,likelihood_only=false;
};
static double number(const std::string& x) {
    size_t p=0; double v=std::stod(x,&p);
    if(p!=x.size() || !std::isfinite(v)) throw std::runtime_error("Invalid numeric argument: "+x);
    return v;
}
static int integer(const std::string& x) {
    double v=number(x);
    if(v<0 || v>100000000 || std::floor(v)!=v) throw std::runtime_error("Invalid integer: "+x);
    return int(v);
}
static std::vector<std::string> split(std::string line) {
    if(!line.empty() && line.back()=='\r') line.pop_back();
    std::vector<std::string> v; std::stringstream ss(line); std::string x;
    while(std::getline(ss,x,'\t')) v.push_back(x);
    return v;
}
static std::ofstream output(const std::string& path) {
    std::ofstream f(path); if(!f) throw std::runtime_error("Cannot write "+path);
    f<<std::setprecision(17); return f;
}
struct Family {std::string id; std::vector<int> counts;};
struct Node {std::string name; double time; int parent=-1,leaf=-1; std::vector<int> children;};
struct Engine {
    std::vector<Node> nodes;
    std::vector<std::string> taxa;
    std::vector<Family> families;
    std::vector<std::pair<std::vector<int>,int>> patterns;
    std::vector<std::vector<std::vector<int>>> subtree_patterns;
    std::vector<int> root_weights;
    Vec prior;
    std::map<double,Vec> matrices;
    int s;
    double lambda,nu,prior_tail=0;
    bool conditioned;
    Engine(const Options& o):s(o.maximum+1),conditioned(!o.unconditioned) {
        std::ifstream f(o.tree); std::string nw; std::getline(f,nw);
        if(!f || nw.find(';')==std::string::npos) throw std::runtime_error("Cannot read Newick tree");
        std::unique_ptr<clade> tree(parse_newick(nw));
        std::map<const clade*,int> indices;
        tree->apply_prefix_order([&](const clade* c){
            int i=nodes.size(); indices[c]=i; Node n;
            n.name=c->is_leaf()?c->get_taxon_name():"Node"+std::to_string(i);
            n.time=c->is_root()?0:c->get_branch_length();
            if(!std::isfinite(n.time)||n.time<0) throw std::runtime_error("Invalid branch length");
            if(!c->is_root()) { n.parent=indices.at(c->get_parent()); nodes[n.parent].children.push_back(i); }
            if(c->is_leaf()) { n.leaf=taxa.size(); taxa.push_back(n.name); }
            nodes.push_back(n);
        });
        if(taxa.size()<2 || std::set<std::string>(taxa.begin(),taxa.end()).size()!=taxa.size())
            throw std::runtime_error("Tree requires at least two distinct tips");
        prior.assign(s,0);
        if(!o.root_file.empty()) {
            std::ifstream r(o.root_file); int k; double p; std::set<int> seen;
            if(!r) throw std::runtime_error("Cannot read root distribution");
            while(r>>k>>p) {
                if(k<0||k>=s||p<0||!std::isfinite(p)||!seen.insert(k).second)
                    throw std::runtime_error("Root prior needs unique counts within --max-count and nonnegative probabilities");
                prior[k]=p;
            }
            if(!r.eof()) throw std::runtime_error("Malformed root distribution; expected count probability pairs, no header");
        } else {
            if(o.root_mean<0) throw std::runtime_error("Specify --root-mean or --root-prior explicitly");
            for(int k=0;k<s;++k) prior[k]=o.root_mean==0?(k==0):std::exp(-o.root_mean+k*std::log(o.root_mean)-std::lgamma(k+1.));
            prior_tail=std::max(0.,1-std::accumulate(prior.begin(),prior.end(),0.));
            if(prior_tail>1e-8) throw std::runtime_error("Root prior tail exceeds 1e-8; increase --max-count");
        }
        double z=std::accumulate(prior.begin(),prior.end(),0.);
        if(z<=0) throw std::runtime_error("Root prior has no positive mass");
        for(double& p:prior) p/=z;
        if(!o.input.empty()) read(o.input);
    }
    void read(const std::string& path) {
        std::ifstream f(path); std::string line; std::getline(f,line); auto head=split(line);
        if(head.size()!=taxa.size()+2) throw std::runtime_error("CAFE table must have Desc, Family ID, and exactly the tree tip columns");
        std::vector<int> cols; std::set<std::string> seen;
        for(size_t j=2;j<head.size();++j) if(!seen.insert(head[j]).second) throw std::runtime_error("Duplicate species column");
        for(auto& taxon:taxa) {
            auto it=std::find(head.begin()+2,head.end(),taxon);
            if(it==head.end()) throw std::runtime_error("Missing taxon: "+taxon);
            cols.push_back(it-head.begin());
        }
        seen.clear(); std::map<std::vector<int>,int> freq;
        while(std::getline(f,line)) {
            if(line.empty()||line=="\r") continue;
            auto fields=split(line);
            if(fields.size()!=head.size()) throw std::runtime_error("Malformed count row");
            Family a; a.id=fields[1];
            if(a.id.empty()||!seen.insert(a.id).second) throw std::runtime_error("Empty or duplicate family ID");
            int sum=0;
            for(int col:cols) {int n=integer(fields[col]); if(n>=s) throw std::runtime_error("Observed count exceeds --max-count: "+a.id); a.counts.push_back(n); sum+=n;}
            if(conditioned && sum==0) throw std::runtime_error("All-zero input family conflicts with observed-family conditioning");
            ++freq[a.counts]; families.push_back(a);
        }
        if(families.empty()) throw std::runtime_error("No families read");
        for(auto& p:freq) patterns.push_back(p);
        // Compile identical subtree observations once. Rates do not enter these keys.
        subtree_patterns.resize(nodes.size());
        std::vector<std::map<std::vector<int>,int>> lookup(nodes.size());
        for(const auto& pattern:patterns) {
            std::vector<int> ids(nodes.size());
            for(int i=int(nodes.size())-1;i>=0;--i) {
                std::vector<int> key;
                if(nodes[i].leaf>=0) key.push_back(pattern.first[nodes[i].leaf]);
                else for(int child:nodes[i].children) key.push_back(ids[child]);
                auto found=lookup[i].find(key);
                if(found==lookup[i].end()) {
                    int id=subtree_patterns[i].size();
                    subtree_patterns[i].push_back(key); lookup[i].emplace(key,id); ids[i]=id;
                    if(i==0) root_weights.push_back(0);
                } else ids[i]=found->second;
            }
            root_weights[ids[0]]+=pattern.second;
        }
    }
    void set_rates(double l,double n) {
        lambda=l; nu=n; matrices.clear();
        for(size_t i=1;i<nodes.size();++i) if(!matrices.count(nodes[i].time))
            matrices.emplace(nodes[i].time,transition(s-1,l,n,nodes[i].time));
    }
    // Scaled pruning: one unit-likelihood vector per node, accumulating log scales.
    double prune(const std::vector<int>& y,std::vector<Vec>* keep=nullptr) const {
        std::vector<Vec> v(nodes.size(),Vec(s,1)); double scale=0;
        for(int i=int(nodes.size())-1;i>=0;--i) {
            const auto& node=nodes[i];
            if(node.leaf>=0) {std::fill(v[i].begin(),v[i].end(),0); v[i][y[node.leaf]]=1;}
            else for(int c:node.children) {
                const auto& a=matrices.at(nodes[c].time);
                for(int j=0;j<s;++j) {
                    double sum=0;
                    if(nodes[c].leaf>=0) sum=a[size_t(j)*s+y[nodes[c].leaf]];
                    else for(int k=0;k<s;++k) sum+=a[size_t(j)*s+k]*v[c][k];
                    v[i][j]*=sum;
                }
                double m=*std::max_element(v[i].begin(),v[i].end());
                if(m<=0) return -INF;
                for(double& x:v[i]) x/=m;
                scale+=std::log(m);
            }
        }
        double z=std::inner_product(prior.begin(),prior.end(),v[0].begin(),0.);
        if(keep) *keep=std::move(v);
        return z>0?std::log(z)+scale:-INF;
    }
    double inclusion_log() const {
        if(!conditioned) return 0;
        double logzero=prune(std::vector<int>(taxa.size(),0));
        if(logzero>=0) return -INF;
        return std::log(-std::expm1(logzero));
    }
    double nll(double l,double n) {
        if(l<0||n<0||!std::isfinite(l)||!std::isfinite(n)) return INF;
        set_rates(l,n); double inc=inclusion_log();
        if(!std::isfinite(inc)) return INF;
        std::vector<std::vector<Vec>> messages(nodes.size());
        std::vector<Vec> scales(nodes.size());
        // A message already includes the transition to the node from its parent.
        // Rebuilt at every (lambda, nu), so cached matrices cannot cross rate values.
        for(int i=int(nodes.size())-1;i>=0;--i) {
            size_t count=subtree_patterns[i].size();
            messages[i].resize(count,Vec(s)); scales[i].resize(count);
            #pragma omp parallel for schedule(static)
            for(size_t p=0;p<count;++p) {
                Vec v(s,1); double logscale=0;
                const auto& key=subtree_patterns[i][p];
                if(nodes[i].leaf>=0) {
                    const auto& a=matrices.at(nodes[i].time);
                    for(int j=0;j<s;++j) messages[i][p][j]=a[size_t(j)*s+key[0]];
                } else {
                    for(size_t c=0;c<nodes[i].children.size();++c) {
                        int child=nodes[i].children[c],id=key[c];
                        logscale+=scales[child][id];
                        for(int j=0;j<s;++j) v[j]*=messages[child][id][j];
                        double m=*std::max_element(v.begin(),v.end());
                        if(m>0) {for(double& x:v) x/=m;logscale+=std::log(m);}
                        else logscale=-INF;
                    }
                    if(i==0) messages[i][p]=std::move(v);
                    else {
                        const auto& a=matrices.at(nodes[i].time);
                        for(int j=0;j<s;++j) {
                            double sum=0;for(int k=0;k<s;++k) sum+=a[size_t(j)*s+k]*v[k];
                            messages[i][p][j]=sum;
                        }
                    }
                }
                double m=*std::max_element(messages[i][p].begin(),messages[i][p].end());
                if(m>0) {for(double& x:messages[i][p]) x/=m;logscale+=std::log(m);}
                else logscale=-INF;
                scales[i][p]=logscale;
            }
        }
        double score=0;
        for(size_t p=0;p<root_weights.size();++p) {
            double z=std::inner_product(prior.begin(),prior.end(),messages[0][p].begin(),0.);
            score-=root_weights[p]*(std::log(z)+scales[0][p]-inc);
        }
        return std::isfinite(score)?score:INF;
    }
    std::vector<Vec> posteriors(const std::vector<int>& y) const {
        std::vector<Vec> inside;
        if(!std::isfinite(prune(y,&inside))) throw std::runtime_error("Cannot reconstruct zero-likelihood family");
        std::vector<Vec> outside(nodes.size(),Vec(s,0)),post=outside; outside[0]=prior;
        for(size_t i=0;i<nodes.size();++i) {
            for(int j=0;j<s;++j) post[i][j]=outside[i][j]*inside[i][j];
            double z=std::accumulate(post[i].begin(),post[i].end(),0.);
            if(z<=0) throw std::runtime_error("Posterior underflow");
            for(double& x:post[i]) x/=z;
            for(int c:nodes[i].children) {
                Vec context=outside[i];
                for(int sibling:nodes[i].children) if(sibling!=c) {
                    const auto& b=matrices.at(nodes[sibling].time);
                    for(int j=0;j<s;++j) {
                        double sum=0; for(int k=0;k<s;++k) sum+=b[size_t(j)*s+k]*inside[sibling][k];
                        context[j]*=sum;
                    }
                    double m=*std::max_element(context.begin(),context.end());
                    if(m>0) for(double& x:context) x/=m;
                }
                const auto& a=matrices.at(nodes[c].time);
                for(int k=0;k<s;++k) for(int j=0;j<s;++j) outside[c][k]+=context[j]*a[size_t(j)*s+k];
                double m=*std::max_element(outside[c].begin(),outside[c].end());
                if(m>0) for(double& x:outside[c]) x/=m;
            }
        }
        return post;
    }
    // Exact, unbounded distribution sampling, independent of transition truncation.
    std::vector<int> sample(std::mt19937& rng) const {
        std::discrete_distribution<int> root(prior.begin(),prior.end());
        for(int attempt=0;attempt<1000000;++attempt) {
            std::vector<int> values(nodes.size()),y(taxa.size()); values[0]=root(rng);
            for(size_t i=1;i<nodes.size();++i) {
                double t=nodes[i].time,x=lambda*t; int parent=values[nodes[i].parent];
                long long value=parent;
                if(lambda>0 && t>0) {
                    int survivors=std::binomial_distribution<int>(parent,1/(1+x))(rng);
                    value=survivors;
                    if(survivors) value+=std::poisson_distribution<int>(std::gamma_distribution<double>(survivors,x)(rng))(rng);
                    if(nu>0) value+=std::poisson_distribution<int>(std::gamma_distribution<double>(nu/lambda,x)(rng))(rng);
                } else if(nu*t>0) value+=std::poisson_distribution<int>(nu*t)(rng);
                if(value>100000000) throw std::runtime_error("Simulation count overflow");
                values[i]=int(value);
            }
            for(size_t i=0;i<nodes.size();++i) if(nodes[i].leaf>=0) y[nodes[i].leaf]=values[i];
            if(!conditioned||std::accumulate(y.begin(),y.end(),0)>0) return y;
        }
        throw std::runtime_error("Observed-family rejection simulation exhausted attempts");
    }
};

struct Scorer:optimizer_scorer {
    Engine& e; double fixed_l,fixed_n; Vec start; int calls=0;
    Scorer(Engine& eng,double l,double n,Vec initial):e(eng),fixed_l(l),fixed_n(n),start(initial) {}
    Vec initial_guesses() override {return start;}
    double calculate_score(const double* x) override {
        int k=0; double l=fixed_l<0?std::exp(x[k++]):fixed_l;
        double n=fixed_n<0?std::exp(x[k++]):fixed_n;
        ++calls;
        if(l>1e4||n>1e4) return INF;
        return e.nll(l,n);
    }
};
struct Fit {double l,n,score; bool converged; int iterations;};
static Fit fit(Engine& e,const Options& o,std::ostream& trace) {
    if(o.lambda>=0 && o.nu>=0) return {o.lambda,o.nu,e.nll(o.lambda,o.nu),true,0};
    Fit best{0,0,INF,false,0};
    // Include exact zero-rate faces, since log parameterization cannot attain zero.
    std::vector<std::pair<double,double>> faces{{o.lambda,o.nu}};
    if(o.lambda<0) faces.push_back({0,o.nu});
    if(o.nu<0) faces.push_back({o.lambda,0});
    if(o.lambda<0 && o.nu<0) faces.push_back({0,0});
    double height=0; for(const auto& node:e.nodes) height=std::max(height,node.time);
    height=std::max(height,1.);
    trace<<"face_lambda\tface_nu\tstart\tlambda\tnu\tnll\titerations\tconverged\n";
    for(auto face:faces) for(int attempt=0;attempt<o.starts;++attempt) {
        Vec initial;
        if(face.first<0) initial.push_back(std::log((.1/height)*std::pow(5.,attempt-1)));
        if(face.second<0) initial.push_back(std::log((.1/height)*std::pow(3.,attempt-1)));
        Fit f{face.first,face.second,INF,true,0};
        if(initial.empty()) {if(attempt) continue; f.score=e.nll(f.l,f.n);f.converged=std::isfinite(f.score);}
        else {
            Scorer scorer(e,face.first,face.second,initial);
            // A near-zero pure-immigration start can underflow for large families.
            // Search for a finite starting likelihood before invoking Nelder-Mead.
            for(int retry=0;retry<12 && !std::isfinite(scorer.calculate_score(initial.data()));++retry)
                for(double& x:initial) x+=1.;
            if(!std::isfinite(scorer.calculate_score(initial.data()))) {
                trace<<face.first<<'\t'<<face.second<<'\t'<<attempt<<"\tNA\tNA\tinf\t0\t0"<<std::endl;
                continue;
            }
            FMinSearch* search=fminsearch_new_with_eq(&scorer,initial.size());
            search->maxiters=o.iterations; search->tolx=1e-5; search->tolf=1e-7;
            fminsearch_min(search,initial.data());
            candidate* c=get_best_result(search); int k=0;
            f.l=face.first<0?std::exp(c->values[k++]):face.first;
            f.n=face.second<0?std::exp(c->values[k++]):face.second;
            f.score=c->score; f.iterations=search->iters; f.converged=std::isfinite(c->score) && search->iters<o.iterations;
            fminsearch_free(search);
        }
        trace<<face.first<<'\t'<<face.second<<'\t'<<attempt<<'\t'<<f.l<<'\t'<<f.n<<'\t'<<f.score<<'\t'<<f.iterations<<'\t'<<f.converged<<std::endl;
        std::cerr<<"BDI fit start: lambda="<<f.l<<" nu="<<f.n<<" nll="<<f.score<<" converged="<<f.converged<<'\n';
        if(f.score<best.score) best=f;
    }
    if(!std::isfinite(best.score)) throw std::runtime_error("No finite likelihood at optimizer starts");
    return best;
}
static void usage() {
    std::cout<<"CAFE5 opt-in innovation model (experimental)\n"
      "cafe5 --innovation -t TREE -i COUNTS --root-mean M -o PREFIX [options]\n"
      "Required root law: --root-mean M (Poisson) OR --root-prior FILE (count probability)\n"
      "--lambda L / --nu N: fix a nonnegative rate; omitted rates are estimated\n"
      "--max-count K (80), --iterations N (400), --starts N (3), --threads N (1)\n"
      "--simulate N --lambda L --nu N: generate observed-family table\n"
      "--likelihood-only: skip family and ancestral output (e.g. likelihood profiles)\n"
      "--seed N; --bootstrap N: experimental fixed-parameter family tail probabilities\n"
      "--unconditioned: disable ascertainment correction AND simulation rejection\n"
      "--no-truncation-check: skip final doubled-state-space likelihood check\n"
      "--matrix-time T --matrix-output FILE --lambda L --nu N: export transition matrix\n"
      "Gamma mixtures, annotation error, branch-specific rates and legacy p-values are unsupported.\n";
}
int run(int argc,char *const argv[]) {
    try {
        Options o;
        for(int i=1;i<argc;++i) {
            std::string a=argv[i]; if(a=="--innovation") continue;
            if(a=="--help"||a=="-h") {usage();return 0;}
            if(a=="--unconditioned") {o.unconditioned=true;continue;}
            if(a=="--no-truncation-check") {o.check=false;continue;}
            if(a=="--likelihood-only") {o.likelihood_only=true;continue;}
            if(i+1>=argc) throw std::runtime_error("Missing value for "+a);
            std::string v=argv[++i];
            if(a=="-t"||a=="--tree") o.tree=v;
            else if(a=="-i"||a=="--infile") o.input=v;
            else if(a=="-o"||a=="--output") o.prefix=v;
            else if(a=="--root-prior") o.root_file=v;
            else if(a=="--root-mean") {o.root_mean=number(v);if(o.root_mean<0) throw std::runtime_error("Negative root mean");}
            else if(a=="--lambda") {o.lambda=number(v);if(o.lambda<0) throw std::runtime_error("Negative lambda");}
            else if(a=="--nu") {o.nu=number(v);if(o.nu<0) throw std::runtime_error("Negative nu");}
            else if(a=="--max-count") o.maximum=integer(v);
            else if(a=="--iterations") o.iterations=integer(v);
            else if(a=="--starts") o.starts=integer(v);
            else if(a=="--simulate") o.simulate=integer(v);
            else if(a=="--bootstrap") o.bootstrap=integer(v);
            else if(a=="--threads") o.threads=integer(v);
            else if(a=="--seed") o.seed=integer(v);
            else if(a=="--matrix-time") o.matrix_time=number(v);
            else if(a=="--matrix-output") o.matrix_file=v;
            else throw std::runtime_error("Unsupported innovation option: "+a);
        }
        if(o.maximum<1||o.maximum>5000||o.starts<1||o.iterations<1||o.threads<1)
            throw std::runtime_error("Invalid computational limit (max-count 1..5000)");
        if(o.root_mean>=0&&!o.root_file.empty()) throw std::runtime_error("Choose one root distribution");
        #ifdef _OPENMP
        omp_set_num_threads(o.threads);
        #endif
        if(!o.matrix_file.empty()) {
            auto a=transition(o.maximum,o.lambda,o.nu,o.matrix_time);auto f=output(o.matrix_file);
            for(int i=0;i<=o.maximum;++i) {for(int j=0;j<=o.maximum;++j) f<<(j?"\t":"")<<a[size_t(i)*(o.maximum+1)+j];f<<'\n';} return 0;
        }
        Engine e(o); std::mt19937 rng(o.seed);
        if(o.simulate) {
            if(o.lambda<0||o.nu<0) throw std::runtime_error("Simulation requires --lambda and --nu");
            e.set_rates(o.lambda,o.nu); if(!std::isfinite(e.inclusion_log())) throw std::runtime_error("Zero inclusion probability");
            auto f=output(o.prefix+"_simulated.tsv");f<<"Desc\tFamily ID";for(auto& t:e.taxa) f<<'\t'<<t;f<<'\n';
            for(int i=0;i<o.simulate;++i) {auto y=e.sample(rng);f<<"BDI\tsim"<<i;for(int n:y) f<<'\t'<<n;f<<'\n';}return 0;
        }
        if(e.families.empty()) throw std::runtime_error("Provide -i COUNTS or --simulate N");
        auto trace=output(o.prefix+"_optimization.tsv"); Fit f=fit(e,o,trace);
        double delta=std::numeric_limits<double>::quiet_NaN();
        if(o.check) {
            Options larger=o; larger.maximum=2*o.maximum;
            Engine check(larger); double large_score=check.nll(f.l,f.n);delta=f.score-large_score;
        }
        e.set_rates(f.l,f.n); double inc=e.inclusion_log();
        auto report=output(o.prefix+"_results.tsv");
        report<<"parameter\tvalue\nmodel\tBDI_equal_birth_death\nlambda\t"<<f.l<<"\nnu\t"<<f.n<<"\nnegative_log_likelihood\t"<<f.score
          <<"\nfamilies\t"<<e.families.size()<<"\nunique_patterns\t"<<e.patterns.size()<<"\nmax_count\t"<<o.maximum
          <<"\nroot_mean\t"<<o.root_mean<<"\nroot_prior_file\t"<<o.root_file<<"\nroot_prior_omitted_mass\t"<<e.prior_tail
          <<"\nconditioned_on_observed\t"<<e.conditioned<<"\nlog_inclusion_probability\t"<<inc
          <<"\noptimizer_converged\t"<<f.converged<<"\ntruncation_nll_difference\t"<<delta
          <<"\ntruncation_pass\t"<<(o.check&&std::isfinite(delta)&&std::abs(delta)<.01)
          <<"\nseed\t"<<o.seed<<"\nsignificance_status\texperimental_only_no_branch_pvalues\n";
        if(o.likelihood_only) {
            if(o.bootstrap) throw std::runtime_error("--likelihood-only cannot be combined with --bootstrap");
            std::cout<<"BDI lambda="<<f.l<<" nu="<<f.n<<" nll="<<f.score<<" truncation_delta="<<delta<<'\n';
            return (f.converged&&(!o.check||(std::isfinite(delta)&&std::abs(delta)<.01)))?0:2;
        }
        auto ancestral=output(o.prefix+"_ancestral.tsv");ancestral<<"Family ID\tNode\tparent\tMAP_count\tposterior_mean\tP_zero\tP_at_cap\n";
        auto likelihood=output(o.prefix+"_families.tsv");likelihood<<"Family ID\tconditional_log_likelihood\n";
        std::map<std::vector<int>,std::vector<Vec>> reconstructed;
        std::vector<double> scores;
        for(auto& family:e.families) {
            double score=e.prune(family.counts)-inc;scores.push_back(score);likelihood<<family.id<<'\t'<<score<<'\n';
            auto it=reconstructed.find(family.counts);
            if(it==reconstructed.end()) it=reconstructed.emplace(family.counts,e.posteriors(family.counts)).first;
            for(size_t i=0;i<e.nodes.size();++i) {
                const auto& p=it->second[i];double mean=0;for(int j=0;j<e.s;++j) mean+=j*p[j];
                ancestral<<family.id<<'\t'<<e.nodes[i].name<<'\t'<<(e.nodes[i].parent<0?"NA":e.nodes[e.nodes[i].parent].name)<<'\t'
                  <<std::distance(p.begin(),std::max_element(p.begin(),p.end()))<<'\t'<<mean<<'\t'<<p[0]<<'\t'<<p.back()<<'\n';
            }
        }
        if(o.bootstrap) {
            Vec null;for(int b=0;b<o.bootstrap;++b) {
                auto y=e.sample(rng);
                if(*std::max_element(y.begin(),y.end())>=e.s) throw std::runtime_error("Bootstrap count exceeds cap: increase --max-count; samples are never silently dropped");
                null.push_back(e.prune(y)-inc);
            }
            std::sort(null.begin(),null.end());auto p=output(o.prefix+"_experimental_pvalues.tsv");
            p<<"Family ID\tfixed_parameter_MC_p\tsimulations\n";
            for(size_t i=0;i<scores.size();++i) p<<e.families[i].id<<'\t'<<(1.+std::distance(null.begin(),std::upper_bound(null.begin(),null.end(),scores[i])))/(o.bootstrap+1.)<<'\t'<<o.bootstrap<<'\n';
        }
        std::cout<<"BDI lambda="<<f.l<<" nu="<<f.n<<" nll="<<f.score<<" truncation_delta="<<delta<<'\n';
        return (f.converged&&(!o.check||(std::isfinite(delta)&&std::abs(delta)<.01)))?0:2;
    } catch(const std::exception& ex) {std::cerr<<"Innovation error: "<<ex.what()<<'\n';return 1;}
}
}
