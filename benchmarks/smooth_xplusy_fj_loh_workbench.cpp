#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <limits>
#include <queue>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>
using ull=unsigned long long;
static constexpr double EPS=1e-12;

static std::vector<int> first_primes(int k){int a[]={2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,83,89}; return std::vector<int>(a,a+k);} 
static std::vector<int> parse_csv(const std::string&s){std::vector<int>v;std::stringstream ss(s);std::string it;while(std::getline(ss,it,','))if(!it.empty())v.push_back(stoi(it));sort(v.begin(),v.end());return v;}
static std::string join_primes(const std::vector<int>&v){std::ostringstream os;os<<"(";for(size_t i=0;i<v.size();i++){if(i)os<<",";os<<v[i];}os<<")";return os.str();}
static std::string idx_to_primes(const std::vector<int>&idx,const std::vector<int>&P){std::ostringstream os;os<<"[";for(size_t i=0;i<idx.size();i++){if(i)os<<",";os<<P[idx[i]];}os<<"]";return os.str();}
static unsigned long long factorial(int k){ unsigned long long r=1; for(int i=2;i<=k;i++) r*=i; return r; }
static double leading_est(const std::vector<int>& primes, ull n){ long double prod=1; for(int p:primes) prod*=log((long double)p); int k=primes.size(); return (double)powl((long double)n*(long double)factorial(k)*prod,1.0L/k); }
static double approx_size(const std::vector<double>& logs,const std::vector<int>& idx,double T){int d=idx.size(); if(!d)return 1; long double prod=1,f=1; for(int i:idx)prod*=logs[i]; for(int j=2;j<=d;j++)f*=j; return (double)(powl(std::max((long double)0.0,(long double)T),d)/(f*prod)+10*powl(std::max((long double)1.0,(long double)T),std::max(0,d-1)));}
struct BuildStats{double sec=0; size_t size=0;};
static std::vector<double> gen_sums_only(const std::vector<double>&logs,const std::vector<int>&idx,double maxT,BuildStats&st){auto t0=std::chrono::high_resolution_clock::now(); int d=idx.size(); std::vector<double> sums; double est=approx_size(logs,idx,maxT); if(est>0 && est<900000000) sums.reserve((size_t)(est*1.02)); sums.push_back(0.0); if(d==0){st={0,1};return sums;} std::vector<size_t> ptr(d,0); std::vector<double>w(d); for(int j=0;j<d;j++)w[j]=logs[idx[j]]; while(true){ double best=1e300; for(int j=0;j<d;j++){ double cand=sums[ptr[j]]+w[j]; if(cand<best) best=cand; } if(best>maxT+EPS) break; sums.push_back(best); for(int j=0;j<d;j++){ while(ptr[j]<sums.size() && sums[ptr[j]]+w[j] <= best+1e-11) ptr[j]++; } } auto t1=std::chrono::high_resolution_clock::now(); st.sec=std::chrono::duration<double>(t1-t0).count(); st.size=sums.size(); return sums;}
static ull pair_count_linear(const std::vector<double>&A,const std::vector<double>&B,double T){long long j=(long long)B.size()-1; ull total=0; for(size_t i=0;i<A.size();i++){double a=A[i]; while(j>=0 && a+B[(size_t)j]>T+EPS)--j; if(j<0)break; total+=(ull)(j+1);} return total;}

struct BlockCounter {
    const std::vector<double>&A; const std::vector<double>&B; size_t leaf; size_t nodes=0, pruned_in=0, pruned_out=0, exact_leaves=0; double T;
    BlockCounter(const std::vector<double>&a,const std::vector<double>&b,size_t leaf_):A(a),B(b),leaf(leaf_){}
    ull exact_slice(size_t a0,size_t a1,size_t b0,size_t b1){
        if(a0>=a1||b0>=b1) return 0; long long j=(long long)b1-1; ull total=0; for(size_t i=a0;i<a1;i++){ double av=A[i]; while(j>=(long long)b0 && av+B[(size_t)j]>T+EPS) --j; if(j<(long long)b0) break; total += (ull)(j-(long long)b0+1); } return total;
    }
    ull rec(size_t a0,size_t a1,size_t b0,size_t b1){
        nodes++;
        if(a0>=a1||b0>=b1) return 0;
        double mn=A[a0]+B[b0]; if(mn>T+EPS){pruned_out++; return 0;}
        double mx=A[a1-1]+B[b1-1]; if(mx<=T+EPS){pruned_in++; return (ull)(a1-a0)*(ull)(b1-b0);} 
        size_t na=a1-a0, nb=b1-b0;
        if(na+nb <= leaf || na<=1 || nb<=1){exact_leaves++; return exact_slice(a0,a1,b0,b1);} 
        if(na>=nb){ size_t mid=a0+na/2; return rec(a0,mid,b0,b1)+rec(mid,a1,b0,b1); }
        else { size_t mid=b0+nb/2; return rec(a0,a1,b0,mid)+rec(a0,a1,mid,b1); }
    }
    ull count(double t){T=t; nodes=pruned_in=pruned_out=exact_leaves=0; return rec(0,A.size(),0,B.size());}
};

struct BlockSel { double sec=0, count_sec=0, band_sec=0, selected=0; int calls=0; ull gap=0, band=0; size_t nodes=0, exact_leaves=0; bool ok=true; };

struct XYData{
  std::vector<int>P; std::vector<double>logs; std::vector<int>Aidx,Bidx; std::vector<double>A,B; BuildStats as,bs; double maxT; double build_sec=0;
  XYData(std::vector<int>p,double T):P(std::move(p)),maxT(T){sort(P.begin(),P.end()); for(int q:P)logs.push_back(log((double)q)); choose_split(); auto t0=std::chrono::high_resolution_clock::now(); A=gen_sums_only(logs,Aidx,maxT,as); B=gen_sums_only(logs,Bidx,maxT,bs); if(A.size()>B.size()){std::swap(A,B);std::swap(Aidx,Bidx);std::swap(as,bs);} auto t1=std::chrono::high_resolution_clock::now(); build_sec=std::chrono::duration<double>(t1-t0).count();}
  double approx(int mask)const{int k=P.size(),d=0;long double prod=1,f=1;for(int i=0;i<k;i++)if(mask>>i&1){d++;prod*=logs[i];}for(int j=2;j<=d;j++)f*=j;return (double)(powl(maxT,d)/(f*prod));}
  void choose_split(){int k=P.size(),full=(1<<k)-1,best=1;double bs=1e300;for(int mask=1;mask<full;mask++){if((mask&1)==0)continue;double a=approx(mask),b=approx(full^mask),score=std::max(a,b); if(score<bs){bs=score;best=mask;}}for(int i=0;i<k;i++)((best>>i)&1?Aidx:Bidx).push_back(i);} 
};

struct AdaptiveSel {double sec=0,count_sec=0,band_sec=0,selected=0;int calls=0;ull gap=0,band=0;};

template<typename CountFn>
AdaptiveSel adaptive_select(const XYData&W, ull n, CountFn&& count_le, ull target_gap=20000){auto t0=std::chrono::high_resolution_clock::now(); double est=leading_est(W.P,n); double lo=std::max(0.0,est*0.70), hi=W.maxT; int calls=0; auto cfn=[&](double t){calls++; return count_le(t);}; ull c_hi=cfn(hi); if(c_hi<n) throw std::runtime_error("maxT too low"); ull c_lo=cfn(lo); while(c_lo>=n){hi=lo;c_hi=c_lo;lo*=.75;c_lo=cfn(lo);} auto tc0=std::chrono::high_resolution_clock::now(); int stagnant=0; for(int iter=0;iter<48;iter++){ull gap=c_hi-c_lo; if(gap<=target_gap)break; long double frac=((long double)n-c_lo)/((long double)c_hi-c_lo); if(!(frac>0&&frac<1))frac=.5; frac=std::min((long double).985,std::max((long double).015,frac)); double t=lo+(double)frac*(hi-lo); if(t<=lo+1e-14*(1+fabs(lo))||t>=hi-1e-14*(1+fabs(hi))||stagnant>=4){t=(lo+hi)/2;stagnant=0;} ull c=cfn(t); ull oldgap=gap; if(c>=n){hi=t;c_hi=c;}else{lo=t;c_lo=c;} ull ng=c_hi-c_lo; if(ng>oldgap*95/100)stagnant++;else stagnant=0;} for(int iter=0;iter<50&&c_hi-c_lo>target_gap;iter++){double mid=(lo+hi)/2; ull c=cfn(mid); if(c>=n){hi=mid;c_hi=c;}else{lo=mid;c_lo=c;}} auto tc1=std::chrono::high_resolution_clock::now(); auto tb0=std::chrono::high_resolution_clock::now(); std::vector<double> vals; vals.reserve((size_t)std::min<ull>(c_hi-c_lo+100,10000000ULL)); long long jhi=(long long)W.B.size()-1, jlo=(long long)W.B.size()-1; for(size_t i=0;i<W.A.size();i++){double a=W.A[i]; if(a>hi+EPS)break; while(jhi>=0 && a+W.B[(size_t)jhi]>hi+EPS)--jhi; if(jhi<0)break; while(jlo>=0 && a+W.B[(size_t)jlo]>lo+EPS)--jlo; long long start=jlo+1, end=jhi; for(long long j=start;j<=end;j++){double s=a+W.B[(size_t)j]; if(s>lo-EPS && s<=hi+EPS) vals.push_back(s);} }
    std::sort(vals.begin(), vals.end()); ull below=pair_count_linear(W.A,W.B,lo); long long off=(long long)n-(long long)below-1; double selected= (off>=0 && (size_t)off<vals.size())? vals[(size_t)off] : hi; auto tb1=std::chrono::high_resolution_clock::now(); auto t1=std::chrono::high_resolution_clock::now(); return {std::chrono::duration<double>(t1-t0).count(), std::chrono::duration<double>(tc1-tc0).count(), std::chrono::duration<double>(tb1-tb0).count(), selected, calls, c_hi-c_lo, (ull)vals.size()};}

struct LOHLayer{size_t s,e; double mn,mx;};
static std::vector<LOHLayer> make_layers_sorted(const std::vector<double>&X, double alpha=1.35){ std::vector<LOHLayer>L; size_t n=X.size(), pos=0; double sz=1.0; while(pos<n){ size_t len=std::max<size_t>(1,(size_t)std::floor(sz)); if(L.size() && len <= L.back().e-L.back().s) len=(L.back().e-L.back().s)+1; size_t e=std::min(n,pos+len); L.push_back({pos,e,X[pos],X[e-1]}); pos=e; sz*=alpha; } return L; }
struct LOHResult{double sec=0; double selected=0; size_t cand=0; size_t layer_pairs=0; bool skipped=false; std::string reason;};
static LOHResult loh_topk_select(const std::vector<double>&A,const std::vector<double>&B,ull k,double alpha=1.35,ull max_cands=30000000ULL){ LOHResult R; auto t0=std::chrono::high_resolution_clock::now(); if(k>max_cands){R.skipped=true;R.reason="top_k_too_large_for_output_style_LOH";return R;} auto LA=make_layers_sorted(A,alpha), LB=make_layers_sorted(B,alpha); struct Tup{double val; int u,v; bool ismax; bool operator>(const Tup&o)const{return val>o.val || (val==o.val && ismax>o.ismax);} }; auto minv=[&](int u,int v){return LA[u].mn+LB[v].mn;}; auto maxv=[&](int u,int v){return LA[u].mx+LB[v].mx;}; std::priority_queue<Tup,std::vector<Tup>,std::greater<Tup>> pq; std::set<std::pair<int,int>> pushed_min, pushed_max; auto push_min=[&](int u,int v){ if(u<0||v<0||u>=(int)LA.size()||v>=(int)LB.size())return; if(pushed_min.insert({u,v}).second) pq.push({minv(u,v),u,v,false});}; auto push_max=[&](int u,int v){ if(pushed_max.insert({u,v}).second) pq.push({maxv(u,v),u,v,true});}; push_min(0,0); ull s=0; std::vector<std::pair<int,int>> q; while(!pq.empty() && s<k){ Tup x=pq.top(); pq.pop(); if(!x.ismax){ push_max(x.u,x.v); push_min(x.u+1,x.v); push_min(x.u,x.v+1); } else { ull prod=(ull)(LA[x.u].e-LA[x.u].s)*(ull)(LB[x.v].e-LB[x.v].s); s += prod; q.push_back({x.u,x.v}); } }
    // Phase 2 from Serang-style LOH selection: include all pending max-tuples,
    // because their layer products can still contain members of the top-k set.
    while(!pq.empty()){ Tup x=pq.top(); pq.pop(); if(x.ismax) q.push_back({x.u,x.v}); }
    std::sort(q.begin(), q.end()); q.erase(std::unique(q.begin(), q.end()), q.end());
    std::vector<double> vals; for(auto [u,v]: q){ size_t na=LA[u].e-LA[u].s, nb=LB[v].e-LB[v].s; if(vals.size()+(ull)na*(ull)nb > max_cands){R.skipped=true;R.reason="candidate_cap_exceeded";return R;} for(size_t i=LA[u].s;i<LA[u].e;i++) for(size_t j=LB[v].s;j<LB[v].e;j++) vals.push_back(A[i]+B[j]); }
    if(vals.size()<k){ R.skipped=true; R.reason="loh_candidate_underflow"; return R; }
    std::nth_element(vals.begin(), vals.begin()+(k-1), vals.end()); R.selected=vals[k-1]; R.cand=vals.size(); R.layer_pairs=q.size(); auto t1=std::chrono::high_resolution_clock::now(); R.sec=std::chrono::duration<double>(t1-t0).count(); return R; }

static ull ceil_quarter(ull x){ return (x+3)/4; }
static ull ma_high_rank(ull x, ull n){ return (n&1) ? ceil_quarter(x+2*n+1) : n+1+ceil_quarter(x); }

struct MARanks {
    long long rank_above_a=0;
    long long rank_at_or_above_b=0;
    std::vector<double> middle;
};

static MARanks ma_ranks_between(const std::vector<double>&X,const std::vector<double>&Y,double a,double b,size_t max_middle){
    ull n=(ull)X.size()-1;
    MARanks r;
    r.rank_at_or_above_b=(long long)(n*n);
    r.middle.reserve((size_t)std::min<ull>(n+1, max_middle));
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

static std::vector<double> ma_half(const std::vector<double>&v){
    std::vector<double> out;
    out.reserve(2+(v.size()-1)/2);
    out.push_back(0.0);
    for(size_t i=1;i<v.size();i+=2) out.push_back(v[i]);
    if(v.size()&1) out.push_back(v.back());
    return out;
}

struct MAPair { double a=0,b=0; };

static MAPair ma_biselect_desc(const std::vector<double>&X,const std::vector<double>&Y,ull k1,ull k2,size_t max_middle){
    ull n=(ull)X.size()-1;
    if(n==0) throw std::runtime_error("empty MA selector input");
    if(n==1){
        double v=X[1]+Y[1];
        return {v,v};
    }
    if(n==2){
        std::vector<double> vals{X[1]+Y[1],X[1]+Y[2],X[2]+Y[1],X[2]+Y[2]};
        std::nth_element(vals.begin(), vals.end()-(long long)k1, vals.end());
        double a=vals[(size_t)(4-k1)];
        std::nth_element(vals.begin(), vals.end()-(long long)k2, vals.end());
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
        std::nth_element(ranks.middle.begin()+1, ranks.middle.end()-r1, ranks.middle.end());
        a=*(ranks.middle.end()-r1);
    }

    double b;
    if(ranks.rank_above_a <= (long long)k2-1) b=coarse.a;
    else if(r2<=0) b=coarse.b;
    else {
        std::nth_element(ranks.middle.begin()+1, ranks.middle.end()-r2, ranks.middle.end());
        b=*(ranks.middle.end()-r2);
    }
    return {a,b};
}

struct MAResult{double sec=0, selected=0; size_t n_square=0, padded_a=0, padded_b=0; bool skipped=false; std::string reason;};

static MAResult ma_select_xplusy_value(const std::vector<double>&A,const std::vector<double>&B,ull k,size_t max_n=30000000ULL,size_t max_middle=80000000ULL){
    MAResult R;
    auto t0=std::chrono::high_resolution_clock::now();
    ull real_pairs=(ull)A.size()*(ull)B.size();
    if(k==0 || k>real_pairs){ R.skipped=true; R.reason="rank_outside_real_product"; return R; }
    size_t n=std::max(A.size(),B.size());
    if(n>max_n){ R.skipped=true; R.reason="square_dimension_cap_exceeded"; return R; }
    const double neg_inf=-std::numeric_limits<double>::infinity();
    std::vector<double>X,Y;
    X.reserve(n+1); Y.reserve(n+1);
    X.push_back(0.0); Y.push_back(0.0);
    for(double x:A) X.push_back(-x);
    for(double y:B) Y.push_back(-y);
    R.padded_a=n-A.size(); R.padded_b=n-B.size();
    while(X.size()<n+1) X.push_back(neg_inf);
    while(Y.size()<n+1) Y.push_back(neg_inf);
    MAPair p=ma_biselect_desc(X,Y,k,k,max_middle);
    R.selected=-p.a;
    R.n_square=n;
    auto t1=std::chrono::high_resolution_clock::now();
    R.sec=std::chrono::duration<double>(t1-t0).count();
    return R;
}

int main(int argc,char**argv){ std::cout.setf(std::ios::fixed); std::cout<<std::setprecision(6); try{
    if(argc>=2 && std::string(argv[1])=="validate-ma"){
        size_t cases=0, failures=0; double max_delta=0.0;
        for(size_t na: {1ul,2ul,3ul,5ul,8ul,13ul,21ul,32ul}){
            for(size_t nb: {1ul,2ul,4ul,7ul,16ul,31ul}){
                std::vector<double>A,B;
                for(size_t i=0;i<na;i++) A.push_back(0.11*(double)i + 0.013*(double)(i*i%7));
                for(size_t j=0;j<nb;j++) B.push_back(0.17*(double)j + 0.019*(double)(j*j%11));
                std::sort(A.begin(),A.end()); std::sort(B.begin(),B.end());
                std::vector<double> vals; vals.reserve(na*nb);
                for(double a:A)for(double b:B)vals.push_back(a+b);
                std::sort(vals.begin(),vals.end());
                for(ull k=1;k<=(ull)vals.size();k++){
                    auto ma=ma_select_xplusy_value(A,B,k,1000,1000000);
                    double delta=ma.skipped?std::numeric_limits<double>::infinity():fabs(ma.selected-vals[(size_t)k-1]);
                    if(delta>max_delta&&std::isfinite(delta))max_delta=delta;
                    if(ma.skipped||delta>1e-10)failures++;
                    cases++;
                }
            }
        }
        std::cout<<"ma_validate cases="<<cases<<" failures="<<failures<<" max_delta="<<std::setprecision(12)<<max_delta<<std::setprecision(6)<<"\n";
        return failures==0?0:1;
    }
    std::vector<std::pair<int,ull>> cases; if(argc>=3){auto P=parse_csv(argv[1]); ull N=std::stoull(argv[2]); cases.push_back({(int)P.size(),N}); double est=leading_est(P,N); double maxT=std::max(1e-9,est*1.02+1e-8); XYData W(P,maxT); auto lin=adaptive_select(W,N,[&](double t){return pair_count_linear(W.A,W.B,t);}); for(size_t leaf: {512ul,2048ul,8192ul}){ BlockCounter bc(W.A,W.B,leaf); auto t0=std::chrono::high_resolution_clock::now(); auto blk=adaptive_select(W,N,[&](double t){return bc.count(t);}); auto t1=std::chrono::high_resolution_clock::now(); std::cout<<"block_rank P="<<join_primes(P)<<" k="<<P.size()<<" N="<<N<<" leaf="<<leaf<<" build="<<W.build_sec<<" A="<<W.A.size()<<" B="<<W.B.size()<<" linear_total="<<lin.sec<<" block_total="<<blk.sec<<" block_calls="<<blk.calls<<" block_log_delta="<<fabs(blk.selected-lin.selected)<<"\n"; }
        auto ma=ma_select_xplusy_value(W.A,W.B,N);
        std::cout<<"ma_select_probe P="<<join_primes(P)<<" k="<<P.size()<<" N="<<N<<" skipped="<<(ma.skipped?"true":"false")<<" sec="<<ma.sec<<" n_square="<<ma.n_square<<" padded_a="<<ma.padded_a<<" padded_b="<<ma.padded_b<<" log="<<std::setprecision(12)<<ma.selected<<std::setprecision(6)<<" log_delta="<<(ma.skipped?0.0:fabs(ma.selected-lin.selected))<<" reason="<<ma.reason<<"\n";
        auto loh=loh_topk_select(W.A,W.B,std::min<ull>(N,1000000ULL),1.35,30000000ULL); std::cout<<"loh_topk_probe P="<<join_primes(P)<<" k="<<P.size()<<" N_probe="<<std::min<ull>(N,1000000ULL)<<" skipped="<<(loh.skipped?"true":"false")<<" sec="<<loh.sec<<" cand="<<loh.cand<<" pairs="<<loh.layer_pairs<<" reason="<<loh.reason<<"\n"; return 0; }

    std::vector<int> ks={4,5,6,8,10,12}; std::vector<ull> Ns={1000000ULL,1000000000ULL,1000000000000ULL}; for(int k:ks){ for(ull N:Ns){ if(k>=12 && N>1000000000ULL) continue; auto P=first_primes(k); double est=leading_est(P,N); double maxT=std::max(1e-9,est*1.02+1e-8); XYData W(P,maxT); auto lin=adaptive_select(W,N,[&](double t){return pair_count_linear(W.A,W.B,t);}); std::cout<<"fj_loh_workbench P="<<join_primes(P)<<" k="<<k<<" N="<<N<<" build="<<W.build_sec<<" A="<<W.A.size()<<" B="<<W.B.size()<<" splitA="<<idx_to_primes(W.Aidx,P)<<" splitB="<<idx_to_primes(W.Bidx,P)<<" linear_total="<<lin.sec<<" linear_count="<<lin.count_sec<<" linear_calls="<<lin.calls<<" linear_band="<<lin.band<<" linear_log="<<std::setprecision(12)<<lin.selected<<std::setprecision(6);
            // block rank test only for not too huge smaller side to keep runtime controlled
            size_t leaf=2048; BlockCounter bc(W.A,W.B,leaf); auto blk=adaptive_select(W,N,[&](double t){return bc.count(t);}); std::cout<<" block_total="<<blk.sec<<" block_count="<<blk.count_sec<<" block_calls="<<blk.calls<<" block_log_delta="<<fabs(blk.selected-lin.selected);
            if(N<=1000000ULL){ auto loh=loh_topk_select(W.A,W.B,N,1.35,30000000ULL); std::cout<<" loh_skipped="<<(loh.skipped?"true":"false")<<" loh_sec="<<loh.sec<<" loh_cand="<<loh.cand<<" loh_pairs="<<loh.layer_pairs; if(!loh.skipped) std::cout<<" loh_delta="<<fabs(loh.selected-lin.selected); else std::cout<<" loh_reason="<<loh.reason; }
            std::cout<<"\n";
        }}
  }catch(std::exception&e){std::cerr<<"error: "<<e.what()<<"\n"; return 1;} }
