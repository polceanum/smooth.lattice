#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include <boost/multiprecision/cpp_int.hpp>

using ull = unsigned long long;
using boost::multiprecision::cpp_int;
static constexpr double EPS = 1e-12;

struct Pack { uint64_t lo=0, hi=0; };
struct Group { std::vector<double> sums; std::vector<Pack> packs; };
struct BuildStats { double sec=0; size_t size=0; };
struct Candidate { double sum; std::vector<int> exps; };

static std::vector<int> parse_csv(const std::string& s) {
    std::vector<int> out; std::stringstream ss(s); std::string item;
    while (std::getline(ss,item,',')) if(!item.empty()) out.push_back(std::stoi(item));
    std::sort(out.begin(), out.end()); return out;
}
static std::string join_vec(const std::vector<int>& v) {
    std::ostringstream os; os<<"[";
    for(size_t i=0;i<v.size();++i){ if(i) os<<","; os<<v[i]; }
    os<<"]"; return os.str();
}
static std::string join_primes(const std::vector<int>& v) {
    std::ostringstream os; os<<"(";
    for(size_t i=0;i<v.size();++i){ if(i) os<<","; os<<v[i]; }
    os<<")"; return os.str();
}
static std::string idx_to_primes(const std::vector<int>& idx, const std::vector<int>& P) {
    std::ostringstream os; os<<"[";
    for(size_t i=0;i<idx.size();++i){ if(i) os<<","; os<<P[idx[i]]; }
    os<<"]"; return os.str();
}
static cpp_int ipow(int base,int exp) {
    cpp_int r=1,b=base;
    while(exp){ if(exp&1) r*=b; exp>>=1; if(exp) b*=b; }
    return r;
}
static cpp_int value_from_exps(const std::vector<int>& primes, const std::vector<int>& exps) {
    cpp_int r=1;
    for(size_t i=0;i<primes.size();++i) if(exps[i]) r*=ipow(primes[i], exps[i]);
    return r;
}
static int digits10(const cpp_int& x) { std::stringstream ss; ss<<x; return (int)ss.str().size(); }
static ull factorial(int k) { ull r=1; for(int i=2;i<=k;i++) r*=i; return r; }
static double leading_est(const std::vector<int>& primes, ull n) {
    double prod=1.0;
    for(int p:primes) prod*=std::log((double)p);
    int k=(int)primes.size();
    return std::pow((double)n*(double)factorial(k)*prod, 1.0/(double)k);
}
static long double factorial_ld(int k) {
    long double r=1.0L;
    for(int i=2;i<=k;i++) r*=(long double)i;
    return r;
}
static std::vector<long double> asymptotic_series_coeffs(const std::vector<int>& primes,int order) {
    std::vector<long double> coeff(order+1,0.0L);
    coeff[0]=1.0L;
    for(int p:primes){
        long double w=std::log((long double)p);
        std::vector<long double> f(order+1,0.0L);
        f[0]=1.0L;
        if(order>=1) f[1]=w/2.0L;
        if(order>=2) f[2]=w*w/12.0L;
        if(order>=4) f[4]=-w*w*w*w/720.0L;
        std::vector<long double> next(order+1,0.0L);
        for(int i=0;i<=order;i++) for(int j=0;j+i<=order;j++) next[i+j]+=coeff[i]*f[j];
        coeff.swap(next);
    }
    long double prod=1.0L;
    for(int p:primes) prod*=std::log((long double)p);
    for(long double& c:coeff) c/=prod;
    return coeff;
}
static long double asymptotic_count_at(const std::vector<long double>& coeff,int k,long double T) {
    long double total=0.0L;
    int order=(int)coeff.size()-1;
    for(int r=0;r<=order && r<=k;r++) total += coeff[r]*std::pow(T,(long double)(k-r))/factorial_ld(k-r);
    return total;
}
static long double asymptotic_count_deriv_at(const std::vector<long double>& coeff,int k,long double T) {
    long double total=0.0L;
    int order=(int)coeff.size()-1;
    for(int r=0;r<=order && r<k;r++) total += coeff[r]*std::pow(T,(long double)(k-r-1))/factorial_ld(k-r-1);
    return total;
}
static double asymptotic_est(const std::vector<int>& primes,ull n) {
    int k=(int)primes.size();
    int order=std::min(k,4);
    auto coeff=asymptotic_series_coeffs(primes,order);
    long double T=leading_est(primes,n);
    long double target=(long double)n;
    for(int iter=0;iter<20;iter++){
        long double f=asymptotic_count_at(coeff,k,T)-target;
        long double d=asymptotic_count_deriv_at(coeff,k,T);
        if(!(d>0.0L) || !std::isfinite((double)d)) break;
        long double step=f/d;
        if(std::fabs(step)>0.25L*T) step=(step>0?0.25L:-0.25L)*T;
        T-=step;
        if(!(T>0.0L) || !std::isfinite((double)T)){ T=leading_est(primes,n); break; }
        if(std::fabs(step)<=1e-13L*(1.0L+T)) break;
    }
    return (double)T;
}
static long double analytic_count_derivative(const std::vector<int>& primes,long double T) {
    int k=(int)primes.size();
    int order=std::min(k,4);
    auto coeff=asymptotic_series_coeffs(primes,order);
    return asymptotic_count_deriv_at(coeff,k,T);
}
static std::string sci(long double x) {
    std::ostringstream os;
    os<<std::scientific<<std::setprecision(12)<<(double)x;
    return os.str();
}
static std::string u128str(__uint128_t x) {
    if(x==0) return "0";
    std::string s;
    while(x){ int d=(int)(x%10); s.push_back('0'+d); x/=10; }
    std::reverse(s.begin(), s.end()); return s;
}
static Pack pack_set(Pack p,int pos,int exp) {
    if(exp<0 || exp>65535) throw std::runtime_error("pack overflow");
    if(pos<4){
        uint64_t mask=0xFFFFULL<<(16*pos);
        p.lo=(p.lo&~mask)|(uint64_t(exp)<<(16*pos));
    } else if(pos<8) {
        int q=pos-4; uint64_t mask=0xFFFFULL<<(16*q);
        p.hi=(p.hi&~mask)|(uint64_t(exp)<<(16*q));
    } else {
        throw std::runtime_error("pack supports <=8 local coords");
    }
    return p;
}
static int pack_get(Pack p,int pos) {
    if(pos<4) return int((p.lo>>(16*pos))&0xFFFFu);
    if(pos<8) return int((p.hi>>(16*(pos-4)))&0xFFFFu);
    throw std::runtime_error("pack supports <=8 local coords");
}
static Pack pack_inc(Pack p,int pos) { return pack_set(p,pos,pack_get(p,pos)+1); }
static bool eligible_parent(Pack pack,int j,int d) {
    for(int h=j+1; h<d; ++h) if(pack_get(pack,h)!=0) return false;
    return true;
}
static double approx_group_size(const std::vector<double>& logs,const std::vector<int>& idx,double T) {
    int d=(int)idx.size();
    if(d==0) return 1.0;
    double prod=1.0;
    for(int i:idx) prod*=logs[i];
    double f=1.0;
    for(int i=2;i<=d;i++) f*=i;
    return std::pow(std::max(0.0,T),d)/(f*prod) + 10.0*std::pow(std::max(1.0,T),std::max(0,d-1));
}
static Group gen_group_pointer_sorted(const std::vector<double>& logs,const std::vector<int>& idx,double maxT,BuildStats& st) {
    auto t0=std::chrono::high_resolution_clock::now();
    Group g;
    double est=approx_group_size(logs,idx,maxT);
    if(est>0 && est<900000000.0){
        g.sums.reserve((size_t)(est*1.03));
        g.packs.reserve((size_t)(est*1.03));
    }
    int d=(int)idx.size();
    g.sums.push_back(0.0); g.packs.push_back(Pack{});
    if(d==0){ st.size=1; st.sec=0; return g; }
    std::vector<size_t> ptr(d,0);
    std::vector<double>w(d);
    for(int j=0;j<d;j++) w[j]=logs[idx[j]];
    while(true){
        double best=1e300; int bestj=-1;
        for(int j=0;j<d;j++){
            while(ptr[j]<g.sums.size() && !eligible_parent(g.packs[ptr[j]],j,d)) ptr[j]++;
            if(ptr[j]>=g.sums.size()) continue;
            double cand=g.sums[ptr[j]]+w[j];
            if(cand<best){ best=cand; bestj=j; }
        }
        if(bestj<0 || best>maxT+EPS) break;
        Pack p=pack_inc(g.packs[ptr[bestj]],bestj);
        g.sums.push_back(best);
        g.packs.push_back(p);
        ptr[bestj]++;
    }
    auto t1=std::chrono::high_resolution_clock::now();
    st.sec=std::chrono::duration<double>(t1-t0).count();
    st.size=g.sums.size();
    return g;
}
static ull pair_count(const std::vector<double>& A,const std::vector<double>& B,double T) {
    long long j=(long long)B.size()-1;
    ull total=0;
    for(size_t i=0;i<A.size();i++){
        double a=A[i];
        while(j>=0 && a+B[(size_t)j]>T+EPS)--j;
        if(j<0)break;
        total+=(ull)(j+1);
    }
    return total;
}

static ull ceil_quarter(ull x){ return (x+3)/4; }
static ull ma_high_rank(ull x, ull n){ return (n&1) ? ceil_quarter(x+2*n+1) : n+1+ceil_quarter(x); }

struct MARanks {
    long long rank_above_a=0;
    long long rank_at_or_above_b=0;
    std::vector<double> middle;
};

static MARanks ma_ranks_between(const std::vector<double>& X,const std::vector<double>& Y,double a,double b,size_t max_middle) {
    ull n=(ull)X.size()-1;
    MARanks r;
    r.rank_at_or_above_b=(long long)(n*n);
    r.middle.reserve((size_t)std::min<ull>(n+1,max_middle));
    r.middle.push_back(0.0);
    long long j=(long long)n;
    for(ull i=1;i<=n;i++){
        while(j>0 && X[(size_t)i]+Y[(size_t)j] <= a) --j;
        r.rank_above_a += j;
        long long jj=j;
        while(jj>0 && X[(size_t)i]+Y[(size_t)jj] < b){
            if(r.middle.size()>=max_middle) throw std::runtime_error("ma_middle_cap_exceeded");
            r.middle.push_back(X[(size_t)i]+Y[(size_t)jj]);
            --jj;
        }
        r.rank_at_or_above_b -= jj;
    }
    return r;
}

static std::vector<double> ma_half(const std::vector<double>& v) {
    std::vector<double> out;
    out.reserve(2+(v.size()-1)/2);
    out.push_back(0.0);
    for(size_t i=1;i<v.size();i+=2) out.push_back(v[i]);
    if(v.size()&1) out.push_back(v.back());
    return out;
}

struct MAPair { double a=0,b=0; };

static MAPair ma_biselect_desc(const std::vector<double>& X,const std::vector<double>& Y,ull k1,ull k2,size_t max_middle) {
    ull n=(ull)X.size()-1;
    if(n==0) throw std::runtime_error("empty MA selector input");
    if(n==1){
        double v=X[1]+Y[1];
        return {v,v};
    }
    if(n==2){
        std::vector<double> vals{X[1]+Y[1],X[1]+Y[2],X[2]+Y[1],X[2]+Y[2]};
        std::nth_element(vals.begin(),vals.end()-(long long)k1,vals.end());
        double a=vals[(size_t)(4-k1)];
        std::nth_element(vals.begin(),vals.end()-(long long)k2,vals.end());
        double b=vals[(size_t)(4-k2)];
        return {a,b};
    }
    ull k1_half=ma_high_rank(k1,n);
    ull k2_half=ceil_quarter(k2);
    MAPair coarse=ma_biselect_desc(ma_half(X),ma_half(Y),k1_half,k2_half,max_middle);
    MARanks ranks=ma_ranks_between(X,Y,coarse.a,coarse.b,max_middle);
    long long nn=(long long)(n*n);
    long long r1=(long long)k1 + ranks.rank_at_or_above_b - nn;
    long long r2=(long long)k2 + ranks.rank_at_or_above_b - nn;

    double a;
    if(ranks.rank_above_a <= (long long)k1-1) a=coarse.a;
    else if(r1<=0) a=coarse.b;
    else {
        std::nth_element(ranks.middle.begin()+1,ranks.middle.end()-r1,ranks.middle.end());
        a=*(ranks.middle.end()-r1);
    }

    double b;
    if(ranks.rank_above_a <= (long long)k2-1) b=coarse.a;
    else if(r2<=0) b=coarse.b;
    else {
        std::nth_element(ranks.middle.begin()+1,ranks.middle.end()-r2,ranks.middle.end());
        b=*(ranks.middle.end()-r2);
    }
    return {a,b};
}

struct MAValueResult {
    double sec=0, selected=0;
    size_t n_square=0, padded_a=0, padded_b=0;
    bool skipped=false;
    std::string reason;
};

static MAValueResult ma_select_xplusy_value(const std::vector<double>& A,const std::vector<double>& B,ull k,size_t max_n,size_t max_middle) {
    MAValueResult r;
    auto t0=std::chrono::high_resolution_clock::now();
    __uint128_t real_pairs=(__uint128_t)A.size()*(__uint128_t)B.size();
    if(k==0 || (__uint128_t)k>real_pairs){ r.skipped=true; r.reason="rank_outside_real_product"; return r; }
    size_t n=std::max(A.size(),B.size());
    if(n>max_n){ r.skipped=true; r.reason="square_dimension_cap_exceeded"; return r; }
    const double neg_inf=-std::numeric_limits<double>::infinity();
    std::vector<double> X,Y;
    X.reserve(n+1); Y.reserve(n+1);
    X.push_back(0.0); Y.push_back(0.0);
    for(double x:A) X.push_back(-x);
    for(double y:B) Y.push_back(-y);
    r.padded_a=n-A.size();
    r.padded_b=n-B.size();
    while(X.size()<n+1) X.push_back(neg_inf);
    while(Y.size()<n+1) Y.push_back(neg_inf);
    MAPair p=ma_biselect_desc(X,Y,k,k,max_middle);
    r.selected=-p.a;
    r.n_square=n;
    auto t1=std::chrono::high_resolution_clock::now();
    r.sec=std::chrono::duration<double>(t1-t0).count();
    return r;
}

struct FullXYRanker {
    std::vector<int> P;
    std::vector<double> logs;
    std::vector<int> Aidx, Bidx;
    Group A, B;
    BuildStats as, bs;
    double maxT=0.0, build_sec=0.0;

    FullXYRanker(std::vector<int> p,double T):P(std::move(p)),maxT(T) {
        std::sort(P.begin(),P.end());
        for(int q:P) logs.push_back(std::log((double)q));
        choose_split();
        auto t0=std::chrono::high_resolution_clock::now();
        A=gen_group_pointer_sorted(logs,Aidx,maxT,as);
        B=gen_group_pointer_sorted(logs,Bidx,maxT,bs);
        if(A.sums.size()>B.sums.size()){
            std::swap(A,B); std::swap(Aidx,Bidx); std::swap(as,bs);
        }
        auto t1=std::chrono::high_resolution_clock::now();
        build_sec=std::chrono::duration<double>(t1-t0).count();
    }
    double approx(int mask)const {
        int k=(int)P.size(),d=0; double prod=1.0,f=1.0;
        for(int i=0;i<k;i++) if(mask>>i&1){ d++; prod*=logs[i]; }
        for(int j=2;j<=d;j++) f*=j;
        return std::pow(maxT,d)/(f*prod);
    }
    void choose_split() {
        int k=(int)P.size(), full=(1<<k)-1, best=1;
        double best_score=1e300;
        for(int mask=1; mask<full; ++mask){
            if((mask&1)==0) continue;
            double a=approx(mask), b=approx(full^mask), score=std::max(a,b);
            if(score<best_score){ best_score=score; best=mask; }
        }
        for(int i=0;i<k;i++) ((best>>i)&1 ? Aidx : Bidx).push_back(i);
    }
    ull count_le(double T) const { return pair_count(A.sums,B.sums,T); }
    std::vector<int> combine(Pack ap, Pack bp) const {
        std::vector<int> e(P.size(),0);
        for(int i=0;i<(int)Aidx.size();++i) e[Aidx[i]]=pack_get(ap,i);
        for(int i=0;i<(int)Bidx.size();++i) e[Bidx[i]]=pack_get(bp,i);
        return e;
    }
    std::vector<Candidate> band(double low,double high) const {
        std::vector<Candidate> out;
        if(high<-EPS) return out;
        long long jhi=(long long)B.sums.size()-1, jlo=(long long)B.sums.size()-1;
        for(size_t ai=0; ai<A.sums.size(); ++ai){
            double a=A.sums[ai];
            if(a>high+EPS) break;
            while(jhi>=0 && a+B.sums[(size_t)jhi]>high+EPS) --jhi;
            if(jhi<0) break;
            while(jlo>=0 && a+B.sums[(size_t)jlo]>low+EPS) --jlo;
            long long start=jlo+1, end=jhi;
            for(long long bj=start; bj<=end; ++bj){
                double s=a+B.sums[(size_t)bj];
                if(s>low-EPS && s<=high+EPS) out.push_back({s,combine(A.packs[ai],B.packs[(size_t)bj])});
            }
        }
        std::sort(out.begin(), out.end(), [](const Candidate& x,const Candidate& y){
            if(x.sum!=y.sum) return x.sum<y.sum;
            return x.exps<y.exps;
        });
        return out;
    }
};

struct Result {
    double sec=0, build=0, count=0, band=0, exact=0;
    std::vector<int> exps;
    int digits=0, calls=0;
    size_t A=0, B=0, bandn=0;
    ull gap=0;
    std::vector<int> Aidx, Bidx;
};

struct CorrectedResult {
    double sec=0, build=0, count=0, band=0, exact=0;
    long double raw_T=0, raw_derivative=0, T=0, derivative=0, correction=0, half_width=0;
    ull center_count=0, below=0, above=0, band_count=0;
    bool target_inside=false, enumerated=false, recovered=false;
    std::vector<int> exps;
    int digits=0;
    size_t A=0, B=0, cands=0;
    std::vector<int> Aidx, Bidx;
};

struct MAFullResult {
    double sec=0, build=0, ma=0, count=0, band=0, exact=0;
    double selected_T=0;
    long double derivative=0, half_width=0;
    ull below=0, above=0, band_count=0;
    bool ma_skipped=false, target_inside=false, enumerated=false, recovered=false;
    std::string ma_reason;
    size_t A=0, B=0, cands=0, n_square=0, padded_a=0, padded_b=0;
    std::vector<int> exps;
    int digits=0, attempts=0;
    std::vector<int> Aidx, Bidx;
};

static CorrectedResult corrected_full_xy(std::vector<int>P, ull n, long double rank_radius, ull max_candidates) {
    auto t0=std::chrono::high_resolution_clock::now();
    CorrectedResult r;
    r.raw_T=(long double)asymptotic_est(P,n);
    r.raw_derivative=analytic_count_derivative(P,r.raw_T);
    if(!(r.raw_derivative>0.0L) || !std::isfinite((double)r.raw_derivative)) throw std::runtime_error("bad analytic derivative");
    long double raw_half_width=std::max(rank_radius/r.raw_derivative,1e-12L*(1.0L+r.raw_T));
    long double initial_slack=std::max(0.01L,10.0L*raw_half_width);
    std::unique_ptr<FullXYRanker> R(new FullXYRanker(P,std::max(1e-9,(double)(r.raw_T+initial_slack))));
    r.build+=R->build_sec;
    auto tc0=std::chrono::high_resolution_clock::now();
    r.center_count=R->count_le((double)r.raw_T);
    auto tc1=std::chrono::high_resolution_clock::now();
    r.count+=std::chrono::duration<double>(tc1-tc0).count();
    r.correction=-((long double)r.center_count-(long double)n)/r.raw_derivative;
    r.T=std::max(0.0L,r.raw_T+r.correction);
    r.derivative=analytic_count_derivative(P,r.T);
    if(!(r.derivative>0.0L) || !std::isfinite((double)r.derivative)) throw std::runtime_error("bad corrected analytic derivative");
    r.half_width=std::max(rank_radius/r.derivative,1e-12L*(1.0L+r.T));
    double L=(double)std::max(0.0L,r.T-r.half_width);
    double H=(double)(r.T+r.half_width);
    if(H+1e-8>R->maxT){
        R.reset(new FullXYRanker(P,std::max(1e-9,H+1e-8)));
        r.build+=R->build_sec;
    }
    auto tc2=std::chrono::high_resolution_clock::now();
    r.below=R->count_le(L);
    r.above=R->count_le(H);
    auto tc3=std::chrono::high_resolution_clock::now();
    r.count+=std::chrono::duration<double>(tc3-tc2).count();
    r.band_count=(r.above>=r.below)?(r.above-r.below):0ULL;
    r.target_inside=(r.below<n && r.above>=n);
    r.A=R->A.sums.size();
    r.B=R->B.sums.size();
    r.Aidx=R->Aidx;
    r.Bidx=R->Bidx;
    if(r.target_inside && r.band_count<=max_candidates){
        auto tb0=std::chrono::high_resolution_clock::now();
        auto cands=R->band(L,H);
        auto tb1=std::chrono::high_resolution_clock::now();
        r.enumerated=true;
        r.cands=cands.size();
        r.band=std::chrono::duration<double>(tb1-tb0).count();
        ull off=n-r.below-1;
        if(off<cands.size()){
            auto te0=std::chrono::high_resolution_clock::now();
            struct EC{ Candidate c; cpp_int val; };
            std::vector<EC> exact;
            exact.reserve(cands.size());
            for(auto& c:cands) exact.push_back({c,value_from_exps(P,c.exps)});
            std::sort(exact.begin(),exact.end(),[](const EC& a,const EC& b){ return a.val<b.val; });
            auto chosen=exact[(size_t)off];
            auto te1=std::chrono::high_resolution_clock::now();
            r.exps=chosen.c.exps;
            r.digits=digits10(chosen.val);
            r.exact=std::chrono::duration<double>(te1-te0).count();
            r.recovered=true;
        }
    }
    auto t1=std::chrono::high_resolution_clock::now();
    r.sec=std::chrono::duration<double>(t1-t0).count();
    return r;
}

static Result nth_full_xy(std::vector<int>P, ull n, ull target_gap=50000) {
    auto t0=std::chrono::high_resolution_clock::now();
    double est=leading_est(P,n), maxT=std::max(1e-9,est*1.02+1e-8);
    std::unique_ptr<FullXYRanker> R(new FullXYRanker(P,maxT));
    int calls=0;
    auto count=[&](double T){ calls++; return R->count_le(T); };
    double lo=std::max(0.0,est*0.70), hi=maxT;
    ull c_hi=count(hi);
    while(c_hi<n){ maxT*=1.25; R.reset(new FullXYRanker(P,maxT)); hi=maxT; c_hi=count(hi); }
    ull c_lo=count(lo);
    while(c_lo>=n){ hi=lo; c_hi=c_lo; lo*=0.75; c_lo=count(lo); }
    auto tc0=std::chrono::high_resolution_clock::now();
    int stagnant=0;
    for(int iter=0;iter<48;iter++){
        ull gap=c_hi-c_lo;
        if(gap<=target_gap) break;
        long double frac=((long double)n-c_lo)/((long double)c_hi-c_lo);
        if(!(frac>0&&frac<1)) frac=.5;
        frac=std::min((long double).985,std::max((long double).015,frac));
        double t=lo+(double)frac*(hi-lo);
        if(t<=lo+1e-14*(1+fabs(lo))||t>=hi-1e-14*(1+fabs(hi))||stagnant>=4){ t=(lo+hi)*0.5; stagnant=0; }
        ull c=count(t), oldgap=gap;
        if(c>=n){ hi=t; c_hi=c; } else { lo=t; c_lo=c; }
        ull ng=c_hi-c_lo;
        if(ng>oldgap*95/100) stagnant++; else stagnant=0;
    }
    for(int iter=0;iter<50&&c_hi-c_lo>target_gap;iter++){
        double mid=(lo+hi)*0.5;
        ull c=count(mid);
        if(c>=n){ hi=mid; c_hi=c; } else { lo=mid; c_lo=c; }
    }
    auto tc1=std::chrono::high_resolution_clock::now();
    auto tb0=std::chrono::high_resolution_clock::now();
    std::vector<Candidate> cands;
    ull below=0; long long off=-1;
    double L=lo,H=hi;
    for(int attempt=0;attempt<8;attempt++){
        below=count(L);
        cands=R->band(L,H);
        if(below<=n-1) off=(long long)n-(long long)below-1;
        if(off>=0 && (size_t)off<cands.size()) break;
        double w=H-L;
        if(w<=0) w=1e-10*(1+fabs(H));
        L=std::max(-1e-9,L-w); H+=w;
    }
    auto tb1=std::chrono::high_resolution_clock::now();
    if(!(off>=0 && (size_t)off<cands.size())) throw std::runtime_error("band failure below="+std::to_string(below)+" cands="+std::to_string(cands.size())+" off="+std::to_string(off));
    auto te0=std::chrono::high_resolution_clock::now();
    struct EC{ Candidate c; cpp_int val; };
    std::vector<EC> exact; exact.reserve(cands.size());
    for(auto& c:cands) exact.push_back({c,value_from_exps(P,c.exps)});
    std::sort(exact.begin(),exact.end(),[](const EC& a,const EC& b){ return a.val<b.val; });
    auto chosen=exact[(size_t)off];
    auto te1=std::chrono::high_resolution_clock::now();
    auto t1=std::chrono::high_resolution_clock::now();
    Result r;
    r.sec=std::chrono::duration<double>(t1-t0).count();
    r.build=R->build_sec;
    r.count=std::chrono::duration<double>(tc1-tc0).count();
    r.band=std::chrono::duration<double>(tb1-tb0).count();
    r.exact=std::chrono::duration<double>(te1-te0).count();
    r.exps=chosen.c.exps;
    r.digits=digits10(chosen.val);
    r.A=R->A.sums.size(); r.B=R->B.sums.size(); r.bandn=cands.size(); r.calls=calls; r.gap=c_hi-c_lo;
    r.Aidx=R->Aidx; r.Bidx=R->Bidx;
    return r;
}

static MAFullResult ma_full_xy(std::vector<int>P, ull n, long double rank_radius, ull max_candidates, size_t max_n, size_t max_middle) {
    auto t0=std::chrono::high_resolution_clock::now();
    double est=asymptotic_est(P,n);
    double maxT=std::max(1e-9,est*1.02+1e-8);
    std::unique_ptr<FullXYRanker> R(new FullXYRanker(P,maxT));
    ull c_hi=R->count_le(maxT);
    while(c_hi<n){
        maxT*=1.25;
        R.reset(new FullXYRanker(P,maxT));
        c_hi=R->count_le(maxT);
    }

    MAFullResult r;
    r.build=R->build_sec;
    r.A=R->A.sums.size();
    r.B=R->B.sums.size();
    r.Aidx=R->Aidx;
    r.Bidx=R->Bidx;

    auto ma=ma_select_xplusy_value(R->A.sums,R->B.sums,n,max_n,max_middle);
    r.ma=ma.sec;
    r.selected_T=ma.selected;
    r.ma_skipped=ma.skipped;
    r.ma_reason=ma.reason;
    r.n_square=ma.n_square;
    r.padded_a=ma.padded_a;
    r.padded_b=ma.padded_b;
    if(ma.skipped){
        auto t1=std::chrono::high_resolution_clock::now();
        r.sec=std::chrono::duration<double>(t1-t0).count();
        return r;
    }

    r.derivative=analytic_count_derivative(P,(long double)r.selected_T);
    if(!(r.derivative>0.0L) || !std::isfinite((double)r.derivative)) throw std::runtime_error("bad MA derivative");
    long double base_width=std::max(rank_radius/r.derivative,1e-12L*(1.0L+(long double)std::fabs(r.selected_T)));
    long double width=base_width;
    for(int attempt=0; attempt<16; ++attempt){
        r.attempts=attempt+1;
        double L=(double)std::max(0.0L,(long double)r.selected_T-width);
        double H=(double)((long double)r.selected_T+width);
        if(H+1e-8>R->maxT){
            R.reset(new FullXYRanker(P,std::max(1e-9,H+1e-8)));
            r.build+=R->build_sec;
            r.A=R->A.sums.size();
            r.B=R->B.sums.size();
            r.Aidx=R->Aidx;
            r.Bidx=R->Bidx;
        }
        auto tc0=std::chrono::high_resolution_clock::now();
        r.below=R->count_le(L);
        r.above=R->count_le(H);
        auto tc1=std::chrono::high_resolution_clock::now();
        r.count+=std::chrono::duration<double>(tc1-tc0).count();
        r.band_count=(r.above>=r.below)?(r.above-r.below):0ULL;
        r.target_inside=(r.below<n && r.above>=n);
        r.half_width=width;
        if(!r.target_inside){
            width*=4.0L;
            continue;
        }
        if(r.band_count>max_candidates){
            break;
        }
        auto tb0=std::chrono::high_resolution_clock::now();
        auto cands=R->band(L,H);
        auto tb1=std::chrono::high_resolution_clock::now();
        r.enumerated=true;
        r.cands=cands.size();
        r.band+=std::chrono::duration<double>(tb1-tb0).count();
        ull off=n-r.below-1;
        if(off<cands.size()){
            auto te0=std::chrono::high_resolution_clock::now();
            struct EC{ Candidate c; cpp_int val; };
            std::vector<EC> exact;
            exact.reserve(cands.size());
            for(auto& c:cands) exact.push_back({c,value_from_exps(P,c.exps)});
            std::sort(exact.begin(),exact.end(),[](const EC& a,const EC& b){ return a.val<b.val; });
            auto chosen=exact[(size_t)off];
            auto te1=std::chrono::high_resolution_clock::now();
            r.exps=chosen.c.exps;
            r.digits=digits10(chosen.val);
            r.exact+=std::chrono::duration<double>(te1-te0).count();
            r.recovered=true;
        }
        break;
    }
    auto t1=std::chrono::high_resolution_clock::now();
    r.sec=std::chrono::duration<double>(t1-t0).count();
    return r;
}

int main(int argc,char**argv) {
    std::cout.setf(std::ios::fixed); std::cout<<std::setprecision(6);
    try {
        if(argc>=2 && std::string(argv[1])=="nth"){
            auto P=parse_csv(argv[2]);
            ull n=std::stoull(argv[3]);
            ull gap=argc>=5?std::stoull(argv[4]):50000ULL;
            auto r=nth_full_xy(P,n,gap);
            std::cout<<"xplusy_full_unrank P="<<join_primes(P)<<" k="<<P.size()<<" N="<<n
                     <<" seconds="<<r.sec<<" build="<<r.build<<" count_phase="<<r.count
                     <<" band_phase="<<r.band<<" exact="<<r.exact<<" calls="<<r.calls
                     <<" rank_gap="<<r.gap<<" A="<<r.A<<" B="<<r.B<<" band="<<r.bandn
                     <<" exps="<<join_vec(r.exps)<<" digits="<<r.digits
                     <<" splitA="<<idx_to_primes(r.Aidx,P)<<" splitB="<<idx_to_primes(r.Bidx,P)<<"\n";
            return 0;
        }
        if(argc>=2 && std::string(argv[1])=="analytic-band-corrected"){
            auto P=parse_csv(argv[2]);
            ull n=std::stoull(argv[3]);
            long double rank_radius=argc>=5?std::stold(argv[4]):1000.0L;
            ull max_candidates=argc>=6?std::stoull(argv[5]):200000ULL;
            auto r=corrected_full_xy(P,n,rank_radius,max_candidates);
            std::cout<<"xplusy_corrected_band P="<<join_primes(P)<<" k="<<P.size()<<" N="<<n
                     <<" seconds="<<r.sec<<" build="<<r.build<<" count_phase="<<r.count
                     <<" band_phase="<<r.band<<" exact="<<r.exact
                     <<" raw_T="<<(double)r.raw_T<<" raw_derivative="<<(double)r.raw_derivative
                     <<" center_count="<<r.center_count
                     <<" center_rank_error="<<(double)((long double)r.center_count-(long double)n)
                     <<" correction="<<sci(r.correction)
                     <<" T="<<(double)r.T<<" derivative="<<(double)r.derivative
                     <<" rank_radius="<<(double)rank_radius<<" half_width="<<sci(r.half_width)
                     <<" below="<<r.below<<" above="<<r.above<<" band_count="<<r.band_count
                     <<" target_inside="<<(r.target_inside?"true":"false")
                     <<" enumerated="<<(r.enumerated?"true":"false")
                     <<" recovered="<<(r.recovered?"true":"false")
                     <<" cands="<<r.cands<<" A="<<r.A<<" B="<<r.B
                     <<" splitA="<<idx_to_primes(r.Aidx,P)<<" splitB="<<idx_to_primes(r.Bidx,P);
            if(r.recovered) std::cout<<" exps="<<join_vec(r.exps)<<" digits="<<r.digits;
            std::cout<<"\n";
            return 0;
        }
        if(argc>=2 && std::string(argv[1])=="ma-full"){
            auto P=parse_csv(argv[2]);
            ull n=std::stoull(argv[3]);
            long double rank_radius=argc>=5?std::stold(argv[4]):1000.0L;
            ull max_candidates=argc>=6?std::stoull(argv[5]):200000ULL;
            size_t max_n=argc>=7?(size_t)std::stoull(argv[6]):30000000ULL;
            size_t max_middle=argc>=8?(size_t)std::stoull(argv[7]):80000000ULL;
            auto r=ma_full_xy(P,n,rank_radius,max_candidates,max_n,max_middle);
            std::cout<<"xplusy_ma_full_unrank P="<<join_primes(P)<<" k="<<P.size()<<" N="<<n
                     <<" seconds="<<r.sec<<" build="<<r.build<<" ma_phase="<<r.ma
                     <<" count_phase="<<r.count<<" band_phase="<<r.band<<" exact="<<r.exact
                     <<" selected_T="<<std::setprecision(12)<<r.selected_T<<std::setprecision(6)
                     <<" derivative="<<(double)r.derivative
                     <<" rank_radius="<<(double)rank_radius<<" half_width="<<sci(r.half_width)
                     <<" below="<<r.below<<" above="<<r.above<<" band_count="<<r.band_count
                     <<" target_inside="<<(r.target_inside?"true":"false")
                     <<" enumerated="<<(r.enumerated?"true":"false")
                     <<" recovered="<<(r.recovered?"true":"false")
                     <<" ma_skipped="<<(r.ma_skipped?"true":"false")
                     <<" ma_reason="<<r.ma_reason
                     <<" attempts="<<r.attempts
                     <<" cands="<<r.cands<<" A="<<r.A<<" B="<<r.B
                     <<" n_square="<<r.n_square<<" padded_a="<<r.padded_a<<" padded_b="<<r.padded_b
                     <<" splitA="<<idx_to_primes(r.Aidx,P)<<" splitB="<<idx_to_primes(r.Bidx,P);
            if(r.recovered) std::cout<<" exps="<<join_vec(r.exps)<<" digits="<<r.digits;
            std::cout<<"\n";
            return r.ma_skipped || !r.recovered ? 1 : 0;
        }
        std::cerr<<"usage: "<<argv[0]<<" nth primes_csv N [target_gap]\n"
                 <<"       "<<argv[0]<<" analytic-band-corrected primes_csv N [rank_radius] [max_candidates]\n"
                 <<"       "<<argv[0]<<" ma-full primes_csv N [rank_radius] [max_candidates] [max_n] [max_middle]\n";
        return 2;
    } catch(const std::exception& e) {
        std::cerr<<"error: "<<e.what()<<"\n";
        return 1;
    }
}
