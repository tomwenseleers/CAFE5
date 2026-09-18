#include "innovation.h"
#include "gamma.h"
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

// One-ancestor linear birth/death PGF: p0 + survival*a*z/(1-q*z).
// Stable on both sides of the critical lambda=mu limit.
struct BDCoefficients {double p0,survival,q,a,log_a;};
static BDCoefficients bd_coefficients(double l,double m,double t) {
    double r=l-m,z=r*t;
    if(r==0) {double x=l*t,a=1/(1+x);return {x*a,a,x*a,a,-std::log1p(x)};}
    if(std::abs(z)<.5) {
        // Preserve log(1-q)/lambda when both per-copy rates are tiny but
        // innovation is not gamma-scaled. Subtracting two logs loses this limit.
        double h=std::expm1(z)/r,d=1+l*h;
        return {m*h/d,std::exp(z)/d,l*h/d,1/d,-std::log1p(l*h)};
    }
    if(z>=0) {
        double h=-std::expm1(-z)/r,d=std::exp(-z)+l*h;
        return {m*h/d,1/d,l*h/d,std::exp(-z)/d,-z-std::log(d)};
    }
    double h=std::expm1(z)/r,d=1+l*h;
    return {m*h/d,std::exp(z)/d,l*h/d,1/d,-std::log1p(l*h)};
}
std::vector<double> transition_asymmetric(int maximum,double l,double m,double nu,double t) {
    if(maximum<0||l<0||m<0||nu<0||t<0||!std::isfinite(l)||!std::isfinite(m)||
       !std::isfinite(nu)||!std::isfinite(t)||!std::isfinite((l+m+nu)*t))
        throw std::runtime_error("Invalid asymmetric BDI transition argument");
    if(l==m)return transition(maximum,l,nu,t);
    int size=maximum+1;Vec matrix(size_t(size)*size,0);
    BDCoefficients c=bd_coefficients(l,m,t);
    double immigration_mean=l==0?nu*(m==0?t:-std::expm1(-m*t)/m):0;
    double lp=l==0?-immigration_mean:(nu==0?0:nu*(c.log_a/l));
    matrix[0]=std::exp(lp);
    for(int j=1;j<size;++j) {
        double numerator=l==0?immigration_mean:nu*(c.q/l)+(j-1)*c.q;
        lp+=(numerator>0?std::log(numerator):-INF)-std::log(double(j));matrix[j]=std::exp(lp);
    }
    for(int i=1;i<size;++i) {
        double sum=0;
        for(int j=0;j<size;++j) {
            if(j)sum=matrix[size_t(i-1)*size+j-1]+c.q*sum;
            matrix[size_t(i)*size+j]=c.p0*matrix[size_t(i-1)*size+j]+c.survival*c.a*sum;
        }
    }
    return matrix;
}

struct Options {
    std::string tree,input,prefix="innovation",root_file,matrix_file,error_file;
    double lambda=-1,nu=-1,root_mean=-1,matrix_time=-1,alpha=-1,epsilon=0;
    double initial_lambda=-1,initial_nu=-1,initial_alpha=-1,initial_epsilon=-1;
    double mu=-1,initial_mu=-1,epsilon_zero=-1,root_zero=.2,root_shape=1;
    bool estimate_mu=false,estimate_epsilon_zero=false,estimate_root_zero=false,estimate_root_shape=false;
    std::string root_family="poisson";
    int categories=1; bool estimate_epsilon=false,epsilon_set=false;
    int maximum=80,iterations=400,starts=3,simulate=0,bootstrap=0,threads=1;
    unsigned seed=20260917;
    bool unconditioned=false,check=true,likelihood_only=false,estimate_root_mean=false,boundary_fits=true;
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
static Vec root_probabilities(int size,const Options& o,double mean,double zero,double shape) {
    Vec p(size,0);
    if(mean<0||!std::isfinite(mean)||zero<0||zero>=1||!std::isfinite(zero)||shape<=0||!std::isfinite(shape))return p;
    if(o.root_family=="poisson") {
        for(int n=0;n<size;++n)p[n]=mean==0?(n==0):std::exp(-mean+n*std::log(mean)-std::lgamma(n+1.));
    } else {
        p[0]=zero;
        for(int n=1;n<size;++n) {
            int k=n-1;double lp;
            if(mean==0)lp=k==0?0:-INF;
            else if(o.root_family=="hurdle-poisson")lp=-mean+k*std::log(mean)-std::lgamma(k+1.);
            else lp=std::lgamma(k+shape)-std::lgamma(shape)-std::lgamma(k+1.)
                    -shape*std::log1p(mean/shape)+k*(std::log(mean)-std::log(shape+mean));
            p[n]=(1-zero)*std::exp(lp);
        }
    }
    return p;
}
// Observation law: P(observed count | true count). At zero, negative counts
// have exactly zero mass; omitted fixed-file rows repeat the preceding row.
struct ErrorLaw {
    double epsilon=0,epsilon_zero=-1;
    double zero_error() const {return epsilon_zero<0?epsilon:epsilon_zero;}
    std::vector<int> deltas{-1,0,1};
    std::map<int,Vec> rows;
    explicit ErrorLaw(const std::string& path="") {
        if(path.empty()) return;
        std::ifstream f(path);if(!f) throw std::runtime_error("Cannot read error model");
        std::string line; bool have_deltas=false;
        while(std::getline(f,line)) {
            if(line.empty()||line=="\r"||line[0]=='#') continue;
            if(line.find("maxcnt:")==0) continue;
            if(line.find("cntdiff:")==0) {
                std::istringstream ss(line.substr(8));int d; deltas.clear();
                while(ss>>d) deltas.push_back(d);
                if(deltas.empty()||!ss.eof()||std::set<int>(deltas.begin(),deltas.end()).size()!=deltas.size())
                    throw std::runtime_error("Invalid error deviations");
                have_deltas=true;continue;
            }
            if(!have_deltas) throw std::runtime_error("Expected cntdiff before error rows");
            std::istringstream ss(line);int count;double v;Vec row;
            if(!(ss>>count)||count<0) throw std::runtime_error("Invalid error row count");
            while(ss>>v) {if(!std::isfinite(v)||v<0) throw std::runtime_error("Invalid error probability");row.push_back(v);}
            if(!ss.eof()||row.size()!=deltas.size()||std::abs(std::accumulate(row.begin(),row.end(),0.)-1)>1e-10||!rows.emplace(count,row).second)
                throw std::runtime_error("Malformed, unnormalized or duplicate error row");
        }
        if(!rows.count(0)) throw std::runtime_error("Error model must explicitly define true count zero");
        int limit=rows.rbegin()->first;
        for(int d:deltas) limit=std::max(limit,-d);
        for(int n=0;n<=limit;++n) {auto row=probabilities(n);for(size_t j=0;j<deltas.size();++j)
            if(n+deltas[j]<0 && row[j]>0) throw std::runtime_error("Error model assigns mass to negative observed counts");}
    }
    Vec probabilities(int truth) const {
        if(rows.empty()) return truth==0?Vec{0,1-zero_error(),zero_error()}:Vec{epsilon,1-2*epsilon,epsilon};
        auto it=rows.upper_bound(truth);--it;return it->second;
    }
    double emission(int observed,int truth) const {
        if(rows.empty()) {
            if(observed==truth)return truth==0?1-zero_error():1-2*epsilon;
            if(observed==truth+1)return truth==0?zero_error():epsilon;
            if(truth>0&&observed==truth-1)return epsilon;
            return 0;
        }
        Vec row=probabilities(truth);
        for(size_t j=0;j<deltas.size();++j) if(observed-truth==deltas[j]) return row[j];
        return 0;
    }
    int sample(int truth,std::mt19937& rng) const {
        if(rows.empty()&&epsilon==0&&zero_error()==0)return truth;
        Vec row=probabilities(truth);int j=std::discrete_distribution<int>(row.begin(),row.end())(rng);
        return truth+deltas[j];
    }
};
static double logsum(const Vec& values) {
    double m=*std::max_element(values.begin(),values.end());
    if(m==-INF) return -INF;
    double sum=0;for(double v:values) sum+=std::exp(v-m);
    return m+std::log(sum);
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
    Vec prior,last_pattern_logs;
    ErrorLaw error;
    std::map<double,Vec> matrices;
    int s;
    double lambda,mu,nu,prior_tail=0;
    bool conditioned;
    Engine(const Options& o):error(o.error_file),s(o.maximum+1),conditioned(!o.unconditioned) {
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
            prior=root_probabilities(s,o,o.root_mean,o.root_zero,o.root_shape);
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
    double leaf_message(const Vec& a,int y,int parent) const {
        if(error.rows.empty()&&error.epsilon==0&&error.zero_error()==0)return a[size_t(parent)*s+y];
        double sum=0;
        for(int d:error.deltas) {int truth=y-d;if(truth>=0&&truth<s)
            sum+=a[size_t(parent)*s+truth]*error.emission(y,truth);}
        return sum;
    }
    void set_rates(double l,double n,double m=-1) {
        lambda=l; mu=m<0?l:m; nu=n; matrices.clear();
        for(size_t i=1;i<nodes.size();++i) if(!matrices.count(nodes[i].time))
            matrices.emplace(nodes[i].time,transition_asymmetric(s-1,l,mu,n,nodes[i].time));
    }
    // Scaled pruning: one unit-likelihood vector per node, accumulating log scales.
    double prune(const std::vector<int>& y,std::vector<Vec>* keep=nullptr) const {
        std::vector<Vec> v(nodes.size(),Vec(s,1)); double scale=0;
        for(int i=int(nodes.size())-1;i>=0;--i) {
            const auto& node=nodes[i];
            if(node.leaf>=0) {std::fill(v[i].begin(),v[i].end(),0); for(int d:error.deltas) {int truth=y[node.leaf]-d;if(truth>=0&&truth<s) v[i][truth]=error.emission(y[node.leaf],truth);}}
            else for(int c:node.children) {
                const auto& a=matrices.at(nodes[c].time);
                for(int j=0;j<s;++j) {
                    double sum=0;
                    if(nodes[c].leaf>=0) sum=leaf_message(a,y[nodes[c].leaf],j);
                    else {
                        #pragma omp simd reduction(+:sum)
                        for(int k=0;k<s;++k) sum+=a[size_t(j)*s+k]*v[c][k];
                    }
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
    double nll(double l,double n,double m=-1) {
        if(l<0||n<0||!std::isfinite(l)||!std::isfinite(n)) return INF;
        set_rates(l,n,m); double inc=inclusion_log();
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
                    for(int j=0;j<s;++j) messages[i][p][j]=leaf_message(a,key[0],j);
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
                            double sum=0;
                            #pragma omp simd reduction(+:sum)
                            for(int k=0;k<s;++k) sum+=a[size_t(j)*s+k]*v[k];
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
        double score=0; last_pattern_logs.resize(root_weights.size());
        for(size_t p=0;p<root_weights.size();++p) {
            double z=std::inner_product(prior.begin(),prior.end(),messages[0][p].begin(),0.);
            last_pattern_logs[p]=std::log(z)+scales[0][p];
            score-=root_weights[p]*(last_pattern_logs[p]-inc);
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
                        double sum=0;
                        #pragma omp simd reduction(+:sum)
                        for(int k=0;k<s;++k) sum+=b[size_t(j)*s+k]*inside[sibling][k];
                        context[j]*=sum;
                    }
                    double m=*std::max_element(context.begin(),context.end());
                    if(m>0) for(double& x:context) x/=m;
                }
                const auto& a=matrices.at(nodes[c].time);
                for(int j=0;j<s;++j) {
                    const double weight=context[j];
                    const double* row=&a[size_t(j)*s];
                    #pragma omp simd
                    for(int k=0;k<s;++k) outside[c][k]+=weight*row[k];
                }
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
                double t=nodes[i].time;int parent=values[nodes[i].parent];
                auto c=bd_coefficients(lambda,mu,t);
                int survivors=std::binomial_distribution<int>(parent,std::min(1.,std::max(0.,c.survival)))(rng);
                long long value=survivors;
                auto poisson_sample=[&](double mean)->int {
                    if(!std::isfinite(mean)||mean>1e8||mean<0)throw std::runtime_error("Simulation intensity overflow");
                    return mean==0?0:std::poisson_distribution<int>(mean)(rng);
                };
                if(lambda>0&&t>0) {
                    double scale=c.q/c.a;
                    if(!std::isfinite(scale))throw std::runtime_error("Simulation geometric scale overflow");
                    if(scale>0&&survivors)value+=poisson_sample(std::gamma_distribution<double>(survivors,scale)(rng));
                    if(scale>0&&nu>0)value+=poisson_sample(std::gamma_distribution<double>(nu/lambda,scale)(rng));
                } else if(nu>0) value+=poisson_sample(nu*(mu==0?t:-std::expm1(-mu*t)/mu));
                if(value>100000000)throw std::runtime_error("Simulation count overflow");
                values[i]=int(value);
            }
            for(size_t i=0;i<nodes.size();++i) if(nodes[i].leaf>=0) y[nodes[i].leaf]=error.sample(values[i],rng);
            if(!conditioned||std::accumulate(y.begin(),y.end(),0)>0) return y;
        }
        throw std::runtime_error("Observed-family rejection simulation exhausted attempts");
    }
};

// Each category is shared across the entire family tree. Selection is applied
// AFTER mixing categories and observation errors, not separately per category.
static Options component_options(Options o) {o.unconditioned=true;return o;}
struct Model:Engine {
    Options options;
    bool observed_condition;
    std::vector<std::unique_ptr<Engine>> extra;
    Vec weights,rates;
    double shape=1,eps=0;
    Model(const Options& o):Engine(component_options(o)),options(o),observed_condition(!o.unconditioned),weights(o.categories,1./o.categories),rates(o.categories,1) {
        for(int k=1;k<o.categories;++k) extra.emplace_back(new Engine(component_options(o)));
        if(observed_condition) for(const auto& f:families)
            if(std::accumulate(f.counts.begin(),f.counts.end(),0)==0) throw std::runtime_error("All-zero observed input family");
    }
    Engine* component(size_t k) {return k?extra[k-1].get():static_cast<Engine*>(this);}
    const Engine* component(size_t k) const {return k?extra[k-1].get():static_cast<const Engine*>(this);}
    bool set_root(double mean,double zero,double root_shape) {
        if(mean<0||!std::isfinite(mean)) return false;
        Vec probabilities=root_probabilities(s,options,mean,zero,root_shape);
        double mass=std::accumulate(probabilities.begin(),probabilities.end(),0.);
        if(!std::isfinite(mass)||mass<=0||1-mass>1e-8) return false;
        for(double& p:probabilities) p/=mass;
        for(size_t k=0;k<rates.size();++k) {component(k)->prior=probabilities;component(k)->prior_tail=std::max(0.,1-mass);}
        return true;
    }
    void set_rates(double l,double n,double a=1,double e=0,double m=-1,double e0=-1) {
        shape=a;eps=e;
        if(rates.size()>1) {
            if(a<.05||a>100) throw std::runtime_error("Gamma shape outside supported interval [0.05,100]");
            get_gamma(weights,rates,a);
        }
        for(size_t k=0;k<rates.size();++k) {
            if(!std::isfinite(rates[k])||rates[k]<=0) throw std::runtime_error("Invalid discrete gamma category");
            component(k)->error.epsilon=e;component(k)->error.epsilon_zero=e0;
            component(k)->set_rates(l*rates[k],n,(m<0?l:m)*rates[k]);
        }
    }
    Vec category_logs(const std::vector<int>& y) const {
        Vec logs(rates.size());for(size_t k=0;k<rates.size();++k) logs[k]=std::log(weights[k])+component(k)->prune(y);
        return logs;
    }
    double prune(const std::vector<int>& y) const {return logsum(category_logs(y));}
    double inclusion_log() const {
        if(!observed_condition) return 0;
        double z=prune(std::vector<int>(taxa.size(),0));
        return z<0?std::log(-std::expm1(z)):-INF;
    }
    double nll(double l,double n,double a=1,double e=0,double m=-1,double e0=-1) {
        if(l<0||n<0||e<0||e>=.5||!std::isfinite(l)||!std::isfinite(n)||!std::isfinite(a)||!std::isfinite(e)||
            (rates.size()>1&&(a<.05||a>100))) return INF;
        set_rates(l,n,a,e,m,e0);
        size_t active=l==0&&(m<0||m==0)?1:rates.size();
        for(size_t k=0;k<active;++k) component(k)->nll(l*rates[k],n,(m<0?l:m)*rates[k]);
        double inc;
        if(active==1) {double z=Engine::prune(std::vector<int>(taxa.size(),0));inc=observed_condition?(z<0?std::log(-std::expm1(z)):-INF):0;}
        else inc=inclusion_log();
        if(!std::isfinite(inc))return INF;
        double score=0;
        for(size_t p=0;p<root_weights.size();++p) {
            Vec logs(active);for(size_t k=0;k<active;++k) logs[k]=(active==1?0:std::log(weights[k]))+component(k)->last_pattern_logs[p];
            score-=root_weights[p]*(logsum(logs)-inc);
        }
        return std::isfinite(score)?score:INF;
    }
    Vec category_posterior(const std::vector<int>& y) const {
        Vec logs=category_logs(y);double total=logsum(logs);
        if(!std::isfinite(total)) throw std::runtime_error("Zero mixture likelihood");
        for(double& v:logs)v=std::exp(v-total);
        return logs;
    }
    std::vector<Vec> posteriors(const std::vector<int>& y) const {
        auto posterior=std::vector<Vec>(nodes.size(),Vec(s,0));Vec w=category_posterior(y);
        for(size_t k=0;k<rates.size();++k) if(w[k]>0) {
            auto post=component(k)->posteriors(y);
            for(size_t i=0;i<nodes.size();++i)for(int j=0;j<s;++j)posterior[i][j]+=w[k]*post[i][j];
        }
        return posterior;
    }
    std::vector<int> sample(std::mt19937& rng) const {
        for(int attempt=0;attempt<1000000;++attempt) {
            size_t k=weights.size()==1?0:std::discrete_distribution<int>(weights.begin(),weights.end())(rng);
            auto y=component(k)->sample(rng);
            if(!observed_condition||std::accumulate(y.begin(),y.end(),0)>0)return y;
        }
        throw std::runtime_error("Observed-mixture simulation exhausted attempts");
    }
};
struct Fit {double l,n,a,e,score;bool converged;int iterations;double root_mean;double mu=-1,e0=-1,root_zero=.2,root_shape=1;};
struct Scorer:optimizer_scorer {
    Model& model;Options options;double fixed_l,fixed_n,fixed_e,fixed_m;Vec start;
    Scorer(Model& m,const Options& o,double l,double n,double e,double mu,Vec v):model(m),options(o),fixed_l(l),fixed_n(n),fixed_e(e),fixed_m(mu),start(v) {}
    Vec initial_guesses() override {return start;}
    Fit decode(const double* x) const {
        int k=0;Fit f;f.l=fixed_l<0?std::exp(x[k++]):fixed_l;f.n=fixed_n<0?std::exp(x[k++]):fixed_n;
        f.mu=options.estimate_mu&&fixed_m<0?std::exp(x[k++]):fixed_m;
        f.a=options.categories==1?1:(options.alpha>=0?options.alpha:((fixed_l==0&&(fixed_m==0||(!options.estimate_mu&&fixed_m<0)))?1:std::exp(x[k++])));
        f.e=fixed_e<0?.5/(1+std::exp(-x[k++])):fixed_e;
        f.root_mean=options.estimate_root_mean?std::exp(x[k++]):options.root_mean;
        f.e0=options.estimate_epsilon_zero?1/(1+std::exp(-x[k++])):options.epsilon_zero;
        f.root_zero=options.estimate_root_zero?1/(1+std::exp(-x[k++])):options.root_zero;
        f.root_shape=options.estimate_root_shape?std::exp(x[k++]):options.root_shape;return f;
    }
    double calculate_score(const double* x) override {
        Fit f=decode(x);if(f.l>1e4||f.n>1e4||f.mu>1e4||!std::isfinite(f.mu)||f.e0>=1)return INF;
        if(options.root_file.empty()&&!model.set_root(f.root_mean,f.root_zero,f.root_shape))return INF;
        return model.nll(f.l,f.n,f.a,f.e,f.mu,f.e0);
    }
};
static Fit fit(Model& e,const Options& o,std::ostream& trace) {
    Fit best;best.score=INF;best.root_mean=o.root_mean;
    std::vector<std::pair<double,double>> faces{{o.lambda,o.nu}};
    if(o.boundary_fits&&o.lambda<0)faces.push_back({0,o.nu});
    if(o.boundary_fits&&o.nu<0)faces.push_back({o.lambda,0});
    if(o.boundary_fits&&o.lambda<0&&o.nu<0)faces.push_back({0,0});
    Vec error_faces{o.estimate_epsilon?-1:o.epsilon};if(o.boundary_fits&&o.estimate_epsilon)error_faces.push_back(0);
    Vec mu_faces{o.estimate_mu?-1:o.mu};if(o.boundary_fits&&o.estimate_mu)mu_faces.push_back(0);
    double height=1;for(const auto& node:e.nodes)height=std::max(height,node.time);
    trace<<"face_lambda\tface_nu\tface_epsilon\tstart\tlambda\tnu\talpha\tepsilon\tnll\titerations\tconverged\troot_mean\tmu\tepsilon_zero\troot_zero\troot_shape\n";
    for(auto face:faces)for(double ef:error_faces)for(double mf:mu_faces)for(int attempt=0;attempt<o.starts;++attempt) {
        Vec initial;
        if(face.first<0)initial.push_back(std::log(attempt==0&&o.initial_lambda>0?o.initial_lambda:(.1/height)*std::pow(5.,attempt-1)));
        if(face.second<0)initial.push_back(std::log(attempt==0&&o.initial_nu>0?o.initial_nu:(.1/height)*std::pow(3.,attempt-1)));
        if(o.estimate_mu&&mf<0)initial.push_back(std::log(o.initial_mu>0?o.initial_mu:(.1/height)*std::pow(5.,attempt-1)));
        if(o.categories>1&&o.alpha<0&&!(face.first==0&&(mf==0||(!o.estimate_mu&&mf<0))))initial.push_back(std::log(attempt==0&&o.initial_alpha>0?o.initial_alpha:std::pow(2.,attempt)));
        if(ef<0){double guess=attempt==0&&o.initial_epsilon>0?o.initial_epsilon:.02*std::pow(2.,attempt);guess=std::min(.2,guess);initial.push_back(std::log(guess/(.5-guess)));}
        if(o.estimate_root_mean)initial.push_back(std::log(std::max(.01,o.root_mean)*std::pow(2.,attempt)));
        if(o.estimate_epsilon_zero){double guess=o.epsilon_zero>0?o.epsilon_zero:.001;initial.push_back(std::log(guess/(1-guess)));}
        if(o.estimate_root_zero)initial.push_back(std::log(o.root_zero/(1-o.root_zero)));
        if(o.estimate_root_shape)initial.push_back(std::log(o.root_shape));
        Scorer scorer(e,o,face.first,face.second,ef,mf,initial);Fit f=scorer.decode(initial.data());f.iterations=0;
        if(initial.empty()){if(attempt)continue;f.score=e.nll(f.l,f.n,f.a,f.e,f.mu,f.e0);f.converged=std::isfinite(f.score);}
        else {
            for(int retry=0;retry<12&&!std::isfinite(scorer.calculate_score(initial.data()));++retry) {
                int k=0;if(face.first<0)initial[k++]+=1;if(face.second<0)initial[k++]+=1;
            }
            if(!std::isfinite(scorer.calculate_score(initial.data())))continue;
            FMinSearch* search=fminsearch_new_with_eq(&scorer,initial.size());
            search->maxiters=o.iterations;search->tolx=1e-5;search->tolf=1e-7;
            fminsearch_min(search,initial.data());candidate* c=get_best_result(search);f=scorer.decode(c->values.data());
            f.score=c->score;f.iterations=search->iters;f.converged=std::isfinite(c->score)&&search->iters<o.iterations;fminsearch_free(search);
        }
        trace<<face.first<<'\t'<<face.second<<'\t'<<ef<<'\t'<<attempt<<'\t'<<f.l<<'\t'<<f.n<<'\t'<<f.a<<'\t'<<f.e<<'\t'<<f.score<<'\t'<<f.iterations<<'\t'<<f.converged<<'\t'<<f.root_mean<<'\t'<<(f.mu<0?f.l:f.mu)<<'\t'<<f.e0<<'\t'<<f.root_zero<<'\t'<<f.root_shape<<std::endl;
        std::cerr<<"BDI fit: lambda="<<f.l<<" nu="<<f.n<<" alpha="<<f.a<<" epsilon="<<f.e<<" nll="<<f.score<<" converged="<<f.converged<<'\n';
        if(f.score<best.score)best=f;
    }
    if(!std::isfinite(best.score))throw std::runtime_error("No finite likelihood at optimizer starts");
    return best;
}
static void usage() {
    std::cout<<"CAFE5 opt-in innovation model (experimental)\n"
      "cafe5 --innovation -t TREE -i COUNTS --root-mean M -o PREFIX [options]\n"
      "Required root law: --root-mean M (Poisson) OR --root-prior FILE (count probability)\n"
      "--lambda L / --nu N: fix a nonnegative rate; omitted rates are estimated\n"
      "--mu M OR --estimate-mu: separate loss rate; default equals lambda\n"
      "--initial-mu M: first loss-rate optimizer start\n"
      "--epsilon-zero E OR --estimate-epsilon-zero: separate zero-to-one error (E supplies an initial value if estimated)\n"
      "--root-family poisson|hurdle-poisson|hurdle-nb; hurdle root-mean is mean excess above one\n"
      "--root-zero P / --root-shape K; --estimate-root-zero / --estimate-root-shape\n"
      "--max-count K (80), --iterations N (400), --starts N (3), --threads N (1)\n"
      "--simulate N --lambda L --nu N: generate observed-family table\n"
      "--likelihood-only: skip family and ancestral output (e.g. likelihood profiles)\n"
      "--skip-boundary-fits: positive-parameter refit only; default also searches exact zero faces\n"
      "--estimate-root-mean: jointly fit Poisson root mean; --root-mean supplies its start\n"
      "--initial-lambda L --initial-nu N --initial-alpha A --initial-epsilon E: first optimizer start only\n"
      "--seed N; --bootstrap N: experimental fixed-parameter family tail probabilities\n"
      "--unconditioned: disable ascertainment correction AND simulation rejection\n"
      "--no-truncation-check: skip final doubled-state-space likelihood check\n"
      "--matrix-time T --matrix-output FILE --lambda L --nu N: export transition matrix\n"
      "--gamma-cats K (1); --alpha A (otherwise estimated when K>1)\n"
      "--epsilon E (fixed) OR --estimate-epsilon OR --error-model FILE\n"
      "Gamma scales duplication/loss only; branch-specific rates and legacy p-values are unsupported.\n";
}
int run(int argc,char *const argv[]) {
    try {
        Options o;
        for(int i=1;i<argc;++i) {
            std::string a=argv[i]; if(a=="--innovation") continue;
            if(a=="--help"||a=="-h") {usage();return 0;}
            if(a=="--unconditioned") {o.unconditioned=true;continue;}
            if(a=="--no-truncation-check") {o.check=false;continue;}
            if(a=="--estimate-epsilon") {o.estimate_epsilon=true;continue;}
            if(a=="--likelihood-only") {o.likelihood_only=true;continue;}
            if(a=="--skip-boundary-fits") {o.boundary_fits=false;continue;}
            if(a=="--estimate-mu") {o.estimate_mu=true;continue;}
            if(a=="--estimate-epsilon-zero") {o.estimate_epsilon_zero=true;continue;}
            if(a=="--estimate-root-zero") {o.estimate_root_zero=true;continue;}
            if(a=="--estimate-root-shape") {o.estimate_root_shape=true;continue;}
            if(a=="--estimate-root-mean") {o.estimate_root_mean=true;continue;}
            if(i+1>=argc) throw std::runtime_error("Missing value for "+a);
            std::string v=argv[++i];
            if(a=="-t"||a=="--tree") o.tree=v;
            else if(a=="-i"||a=="--infile") o.input=v;
            else if(a=="-o"||a=="--output") o.prefix=v;
            else if(a=="--root-prior") o.root_file=v;
            else if(a=="--root-mean") {o.root_mean=number(v);if(o.root_mean<0) throw std::runtime_error("Negative root mean");}
            else if(a=="--lambda") {o.lambda=number(v);if(o.lambda<0) throw std::runtime_error("Negative lambda");}
            else if(a=="--nu") {o.nu=number(v);if(o.nu<0) throw std::runtime_error("Negative nu");}
            else if(a=="--mu") {o.mu=number(v);if(o.mu<0)throw std::runtime_error("Negative mu");}
            else if(a=="--initial-mu") o.initial_mu=number(v);
            else if(a=="--epsilon-zero") o.epsilon_zero=number(v);
            else if(a=="--root-family") o.root_family=v;
            else if(a=="--root-zero") o.root_zero=number(v);
            else if(a=="--root-shape") o.root_shape=number(v);
            else if(a=="--initial-lambda") o.initial_lambda=number(v);
            else if(a=="--initial-nu") o.initial_nu=number(v);
            else if(a=="--initial-alpha") o.initial_alpha=number(v);
            else if(a=="--initial-epsilon") o.initial_epsilon=number(v);
            else if(a=="--gamma-cats") o.categories=integer(v);
            else if(a=="--alpha") o.alpha=number(v);
            else if(a=="--epsilon") {o.epsilon=number(v);o.epsilon_set=true;}
            else if(a=="--error-model") o.error_file=v;
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
        if(o.root_family!="poisson"&&o.root_family!="hurdle-poisson"&&o.root_family!="hurdle-nb")throw std::runtime_error("Unknown root family");
        if(o.root_zero<0||o.root_zero>=1||o.root_shape<=0||(o.epsilon_zero<0&&o.epsilon_zero!=-1)||o.epsilon_zero>=1)throw std::runtime_error("Invalid root/error parameter");
        if(o.estimate_mu&&o.mu>=0)throw std::runtime_error("Choose fixed or estimated mu");
        if((o.estimate_root_zero||o.estimate_root_shape||o.root_family!="poisson")&&!o.root_file.empty())throw std::runtime_error("Parametric root options conflict with root file");
        if(o.estimate_root_zero&&(o.root_family=="poisson"||o.root_zero<=0))throw std::runtime_error("Estimated root zero needs a hurdle root and an interior starting value");
        if(o.estimate_root_shape&&o.root_family!="hurdle-nb")throw std::runtime_error("Root shape estimation requires hurdle-nb");
        if(!o.error_file.empty()&&(o.epsilon_zero>=0||o.estimate_epsilon_zero))throw std::runtime_error("Zero-error options conflict with error file");
        if(o.maximum<1||o.maximum>5000||o.starts<1||o.iterations<1||o.threads<1)
            throw std::runtime_error("Invalid computational limit (max-count 1..5000)");
        if(o.estimate_root_mean&&(o.root_mean<=0||!o.root_file.empty()||o.simulate)) throw std::runtime_error("Estimating root mean requires positive --root-mean initial value, no root file, and observed data");
        if(o.root_mean>=0&&!o.root_file.empty()) throw std::runtime_error("Choose one root distribution");
        #ifdef _OPENMP
        omp_set_num_threads(o.threads);
        #endif
        if(!o.matrix_file.empty()) {
            auto a=transition_asymmetric(o.maximum,o.lambda,o.mu<0?o.lambda:o.mu,o.nu,o.matrix_time);auto f=output(o.matrix_file);
            for(int i=0;i<=o.maximum;++i) {for(int j=0;j<=o.maximum;++j) f<<(j?"\t":"")<<a[size_t(i)*(o.maximum+1)+j];f<<'\n';} return 0;
        }
        if(o.categories<1||o.categories>32||o.epsilon<0||o.epsilon>=.5||
            (o.alpha!=-1&&(o.alpha<.05||o.alpha>100))||(o.categories==1&&o.alpha!=-1))
            throw std::runtime_error("Invalid gamma/error settings");
        if(int(o.epsilon_set)+int(o.estimate_epsilon)+int(!o.error_file.empty())>1)
            throw std::runtime_error("Choose fixed epsilon, estimated epsilon OR an error-model file");
        Model e(o); std::mt19937 rng(o.seed);
        if(o.simulate) {
            if(o.lambda<0||o.nu<0) throw std::runtime_error("Simulation requires --lambda and --nu");
            if(o.estimate_mu||o.estimate_epsilon_zero||o.estimate_root_zero||o.estimate_root_shape||o.estimate_epsilon||(o.categories>1&&o.alpha<0))throw std::runtime_error("Simulation requires fixed alpha and epsilon");
            e.set_rates(o.lambda,o.nu,o.categories>1?o.alpha:1,o.epsilon,o.mu,o.epsilon_zero); if(!std::isfinite(e.inclusion_log())) throw std::runtime_error("Zero inclusion probability");
            auto f=output(o.prefix+"_simulated.tsv");f<<"Desc\tFamily ID";for(auto& t:e.taxa) f<<'\t'<<t;f<<'\n';
            for(int i=0;i<o.simulate;++i) {auto y=e.sample(rng);f<<"BDI\tsim"<<i;for(int n:y) f<<'\t'<<n;f<<'\n';}return 0;
        }
        if(e.families.empty()) throw std::runtime_error("Provide -i COUNTS or --simulate N");
        auto trace=output(o.prefix+"_optimization.tsv"); Fit f=fit(e,o,trace);
        double delta=std::numeric_limits<double>::quiet_NaN();
        if(o.check) {
            Options larger=o; larger.maximum=2*o.maximum;
            Model check(larger); if(o.root_file.empty()&&!check.set_root(f.root_mean,f.root_zero,f.root_shape))throw std::runtime_error("Root prior truncation failed"); double large_score=check.nll(f.l,f.n,f.a,f.e,f.mu,f.e0);delta=f.score-large_score;
        }
        if(o.root_file.empty()&&!e.set_root(f.root_mean,f.root_zero,f.root_shape))throw std::runtime_error("Root prior truncation failed");
        e.set_rates(f.l,f.n,f.a,f.e,f.mu,f.e0); double inc=e.inclusion_log();
        auto report=output(o.prefix+"_results.tsv");
        report<<"parameter\tvalue\nmodel\t"<<((o.mu>=0||o.estimate_mu)?"BDI_separate_birth_death":"BDI_equal_birth_death")<<"\nlambda\t"<<f.l<<"\nnu\t"<<f.n<<"\nnegative_log_likelihood\t"<<f.score
          <<"\nfamilies\t"<<e.families.size()<<"\nunique_patterns\t"<<e.patterns.size()<<"\nmax_count\t"<<o.maximum
          <<"\ngamma_categories\t"<<o.categories<<"\nalpha\t"<<f.a<<"\nepsilon\t"<<f.e<<"\nerror_model_file\t"<<o.error_file<<"\nalpha_at_bound\t"<<(o.categories>1&&(f.a<.0501||f.a>99.99))
          <<"\nalpha_identifiable\t"<<(o.categories>1&&(f.l>0||f.mu>0))
          <<"\nmu\t"<<(f.mu<0?f.l:f.mu)<<"\nmu_estimated\t"<<o.estimate_mu<<"\nepsilon_zero\t"<<(f.e0<0?f.e:f.e0)<<"\nepsilon_zero_separate\t"<<(f.e0>=0)<<"\nepsilon_zero_estimated\t"<<o.estimate_epsilon_zero<<"\nroot_zero_estimated\t"<<o.estimate_root_zero<<"\nroot_shape_estimated\t"<<o.estimate_root_shape<<"\nP_root_zero\t"<<e.prior[0]<<"\nroot_family\t"<<o.root_family<<"\nroot_zero\t"<<f.root_zero<<"\nroot_shape\t"<<f.root_shape
          <<"\nroot_mean\t"<<f.root_mean<<"\nroot_mean_estimated\t"<<o.estimate_root_mean<<"\nroot_prior_file\t"<<o.root_file<<"\nroot_prior_omitted_mass\t"<<e.prior_tail
          <<"\nconditioned_on_observed\t"<<e.observed_condition<<"\nlog_inclusion_probability\t"<<inc
          <<"\nboundary_fits\t"<<o.boundary_fits<<"\noptimizer_converged\t"<<f.converged<<"\ntruncation_nll_difference\t"<<delta
          <<"\ntruncation_pass\t"<<(o.check&&std::isfinite(delta)&&std::abs(delta)<.01)
          <<"\nseed\t"<<o.seed<<"\nsignificance_status\tuse_branch_bootstrap_driver_for_calibrated_tail_tests\n";
        if(o.likelihood_only) {
            if(o.bootstrap) throw std::runtime_error("--likelihood-only cannot be combined with --bootstrap");
            std::cout<<"BDI lambda="<<f.l<<" nu="<<f.n<<" nll="<<f.score<<" truncation_delta="<<delta<<'\n';
            return (f.converged&&(!o.check||(std::isfinite(delta)&&std::abs(delta)<.01)))?0:2;
        }
        auto categories=output(o.prefix+"_categories.tsv");categories<<"category\tweight\tlambda_multiplier\n";
        for(size_t k=0;k<e.rates.size();++k)categories<<k<<'\t'<<e.weights[k]<<'\t'<<e.rates[k]<<'\n';
        auto branch=output(o.prefix+"_branch_statistics.tsv");branch<<"Family ID\tNode\tparent\tposterior_mean_change\tparent_MAP_count\tchild_MAP_count\tMAP_change\ttransition_tail_score\n";
        auto ancestral=output(o.prefix+"_ancestral.tsv");ancestral<<"Family ID\tNode\tparent\tMAP_count\tposterior_mean\tP_zero\tP_at_cap\n";
        auto likelihood=output(o.prefix+"_families.tsv");likelihood<<"Family ID\tconditional_log_likelihood";for(size_t k=0;k<e.rates.size();++k)likelihood<<"\tP_category_"<<k;likelihood<<'\n';
        std::map<std::vector<int>,std::vector<Vec>> reconstructed;
        // Compute each distinct posterior once in parallel; output stays deterministic.
        std::vector<std::vector<Vec>> pattern_posteriors(e.patterns.size());
        std::vector<std::string> reconstruction_errors(e.patterns.size());
        std::vector<Vec> pattern_categories(e.patterns.size());
        Vec pattern_scores(e.patterns.size());
        #pragma omp parallel for schedule(dynamic)
        for(size_t k=0;k<e.patterns.size();++k) {
            try {
                pattern_posteriors[k]=e.posteriors(e.patterns[k].first);
                pattern_scores[k]=e.prune(e.patterns[k].first)-inc;
                pattern_categories[k]=e.category_posterior(e.patterns[k].first);
            }
            catch(const std::exception& ex) {reconstruction_errors[k]=ex.what();}
        }
        std::map<std::vector<int>,size_t> pattern_index;
        for(size_t k=0;k<e.patterns.size();++k) {
            pattern_index.emplace(e.patterns[k].first,k);
            if(!reconstruction_errors[k].empty()) throw std::runtime_error(reconstruction_errors[k]);
            reconstructed.emplace(e.patterns[k].first,std::move(pattern_posteriors[k]));
        }
        std::vector<double> scores;
        for(auto& family:e.families) {
            size_t pattern_id=pattern_index.at(family.counts);
            double score=pattern_scores[pattern_id];scores.push_back(score);likelihood<<family.id<<'\t'<<score;for(double w:pattern_categories[pattern_id])likelihood<<'\t'<<w;likelihood<<'\n';
            auto it=reconstructed.find(family.counts);
            if(it==reconstructed.end()) it=reconstructed.emplace(family.counts,e.posteriors(family.counts)).first;
            Vec means(e.nodes.size());for(size_t i=0;i<e.nodes.size();++i)for(int j=0;j<e.s;++j)means[i]+=j*it->second[i][j];
            std::vector<int> modes(e.nodes.size());
            for(size_t i=0;i<e.nodes.size();++i) modes[i]=std::distance(it->second[i].begin(),std::max_element(it->second[i].begin(),it->second[i].end()));
            for(size_t i=1;i<e.nodes.size();++i) {
                int parent=e.nodes[i].parent,from=modes[parent],to=modes[i];
                double lower=0,upper=0;
                for(size_t k=0;k<e.rates.size();++k) {
                    const auto& matrix=e.component(k)->matrices.at(e.nodes[i].time);
                    double before=0;for(int j=0;j<to;++j)before+=matrix[size_t(from)*e.s+j];
                    const double weight=pattern_categories[pattern_id][k];
                    lower+=weight*(before+matrix[size_t(from)*e.s+to]);
                    // The upper tail includes mass beyond the numerical cap.
                    upper+=weight*std::max(0.,1-before);
                }
                double tail=std::min(1.,2*std::min(lower,upper));
                double surprise=-std::log(std::max(std::numeric_limits<double>::min(),tail));
                branch<<family.id<<'\t'<<e.nodes[i].name<<'\t'<<e.nodes[parent].name<<'\t'<<means[i]-means[parent]
                    <<'\t'<<from<<'\t'<<to<<'\t'<<to-from<<'\t'<<surprise<<'\n';
            }
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
