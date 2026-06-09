
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <iomanip>
#include <limits>
#include <iostream>
#include <memory>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>
#include <boost/multiprecision/cpp_int.hpp>
using ull = unsigned long long;
using u128 = __uint128_t;
using boost::multiprecision::cpp_int;
static constexpr double EPS = 1e-12;

static std::vector<int> parse_csv(const std::string& s){ std::vector<int> out; std::stringstream ss(s); std::string item; while(std::getline(ss,item,',')) if(!item.empty()) out.push_back(std::stoi(item)); std::sort(out.begin(), out.end()); return out; }
static std::string join_vec(const std::vector<int>& v){ std::ostringstream os; os<<"["; for(size_t i=0;i<v.size();++i){ if(i) os<<","; os<<v[i]; } os<<"]"; return os.str(); }
static std::string join_primes(const std::vector<int>& v){ std::ostringstream os; os<<"("; for(size_t i=0;i<v.size();++i){ if(i) os<<","; os<<v[i]; } os<<")"; return os.str(); }
static std::string idx_to_primes(const std::vector<int>& idx,const std::vector<int>& P){ std::ostringstream os; os<<"["; for(size_t i=0;i<idx.size();++i){ if(i) os<<","; os<<P[idx[i]]; } os<<"]"; return os.str(); }
static std::string u128str(u128 x){ if(x==0) return "0"; std::string s; while(x){ int d=(int)(x%10); s.push_back('0'+d); x/=10; } std::reverse(s.begin(),s.end()); return s; }
static cpp_int ipow(int base,int exp){ cpp_int r=1,b=base; while(exp){ if(exp&1) r*=b; exp>>=1; if(exp) b*=b; } return r; }
static cpp_int value_from_exps(const std::vector<int>& primes,const std::vector<int>& exps){ cpp_int r=1; for(size_t i=0;i<primes.size();++i) if(exps[i]) r*=ipow(primes[i],exps[i]); return r; }
static int digits10(const cpp_int& x){ std::stringstream ss; ss<<x; return (int)ss.str().size(); }
static unsigned long long factorial(int k){ unsigned long long r=1; for(int i=2;i<=k;i++) r*=i; return r; }
static double leading_est(const std::vector<int>& primes, ull n){ double prod=1.0; for(int p:primes) prod*=std::log((double)p); int k=(int)primes.size(); return std::pow((double)n*(double)factorial(k)*prod,1.0/(double)k); }
static long double factorial_ld(int k){ long double r=1.0L; for(int i=2;i<=k;i++) r*=(long double)i; return r; }
static std::vector<long double> asymptotic_series_coeffs(const std::vector<int>& primes,int order){
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
static long double asymptotic_count_at(const std::vector<long double>& coeff,int k,long double T){
    long double total=0.0L;
    int order=(int)coeff.size()-1;
    for(int r=0;r<=order && r<=k;r++) total += coeff[r]*std::pow(T,(long double)(k-r))/factorial_ld(k-r);
    return total;
}
static long double asymptotic_count_deriv_at(const std::vector<long double>& coeff,int k,long double T){
    long double total=0.0L;
    int order=(int)coeff.size()-1;
    for(int r=0;r<=order && r<k;r++) total += coeff[r]*std::pow(T,(long double)(k-r-1))/factorial_ld(k-r-1);
    return total;
}
static double asymptotic_est(const std::vector<int>& primes,ull n){
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
        if(std::fabsl(step)>0.25L*T) step=(step>0?0.25L:-0.25L)*T;
        T-=step;
        if(!(T>0.0L) || !std::isfinite((double)T)){ T=leading_est(primes,n); break; }
        if(std::fabsl(step)<=1e-13L*(1.0L+T)) break;
    }
    return (double)T;
}

struct Pack{ uint64_t lo=0, hi=0; };
static Pack pack_set(Pack p,int pos,int exp){ if(exp<0||exp>65535) throw std::runtime_error("pack overflow"); if(pos<4){ uint64_t mask=0xFFFFULL<<(16*pos); p.lo=(p.lo&~mask)|(uint64_t(exp)<<(16*pos)); } else if(pos<8){ int q=pos-4; uint64_t mask=0xFFFFULL<<(16*q); p.hi=(p.hi&~mask)|(uint64_t(exp)<<(16*q)); } else throw std::runtime_error("pack supports <=8 local coords"); return p; }
static int pack_get(Pack p,int pos){ if(pos<4) return int((p.lo>>(16*pos))&0xFFFFu); if(pos<8) return int((p.hi>>(16*(pos-4)))&0xFFFFu); throw std::runtime_error("pack supports <=8 local coords"); }
static Pack pack_inc(Pack p,int pos){ return pack_set(p,pos,pack_get(p,pos)+1); }

struct Group{ std::vector<double> sums; std::vector<Pack> packs; };
struct BuildStats{ double sec=0; size_t size=0; };
static double approx_group_size(const std::vector<double>& logs,const std::vector<int>& idx,double T){ int d=(int)idx.size(); if(d==0) return 1.0; double prod=1.0; for(int i:idx) prod*=logs[i]; double f=1.0; for(int i=2;i<=d;i++) f*=i; return std::pow(std::max(0.0,T),d)/(f*prod)+10.0*std::pow(std::max(1.0,T),std::max(0,d-1)); }
static bool eligible_parent(Pack pack,int j,int d){ for(int h=j+1; h<d; ++h) if(pack_get(pack,h)!=0) return false; return true; }
static Group gen_group_pointer_sorted(const std::vector<double>& logs,const std::vector<int>& idx,double maxT,BuildStats& st){
    auto t0=std::chrono::high_resolution_clock::now(); Group g; double est=approx_group_size(logs,idx,maxT); if(est>0 && est<300000000.0){ g.sums.reserve((size_t)(est*1.03)); g.packs.reserve((size_t)(est*1.03)); }
    int d=(int)idx.size(); g.sums.push_back(0.0); g.packs.push_back(Pack{}); if(d==0){ st.size=1; st.sec=0; return g; }
    std::vector<size_t> ptr(d,0); std::vector<double>w(d); for(int j=0;j<d;j++) w[j]=logs[idx[j]];
    while(true){ double best=1e300; int bestj=-1; for(int j=0;j<d;j++){ while(ptr[j]<g.sums.size() && !eligible_parent(g.packs[ptr[j]],j,d)) ptr[j]++; if(ptr[j]>=g.sums.size()) continue; double cand=g.sums[ptr[j]]+w[j]; if(cand<best){ best=cand; bestj=j; } } if(bestj<0 || best>maxT+EPS) break; Pack p=pack_inc(g.packs[ptr[bestj]],bestj); g.sums.push_back(best); g.packs.push_back(p); ptr[bestj]++; }
    auto t1=std::chrono::high_resolution_clock::now(); st.sec=std::chrono::duration<double>(t1-t0).count(); st.size=g.sums.size(); return g;
}
struct Fenwick{ int n; std::vector<int> bit; Fenwick(int n_=0):n(n_),bit(n_+1,0){} void reset(int n_){n=n_;bit.assign(n+1,0);} void add(int i,int v){ for(++i;i<=n;i+=i&-i) bit[i]+=v; } int sum_prefix(int i) const{ int r=0; for(;i>0;i-=i&-i) r+=bit[i]; return r; } int sum_le_index(int idx) const{ if(idx<0) return 0; if(idx>=n) idx=n-1; return sum_prefix(idx+1); } };

struct CompressedLayer{
    const Group* base=nullptr; double w=0.0; std::vector<ull> qval, prefix_q; std::vector<double> resid,resid_sorted; std::vector<int> resid_rank; std::vector<size_t> resid_order;
    CompressedLayer(){} CompressedLayer(const Group& g,double outer_w){ init(g,outer_w); }
    void init(const Group& g,double outer_w){ base=&g; w=outer_w; size_t n=g.sums.size(); qval.resize(n); prefix_q.assign(n+1,0); resid.resize(n); resid_sorted.resize(n); resid_rank.resize(n); resid_order.resize(n); for(size_t i=0;i<n;i++){ double x=g.sums[i]; ull q=(ull)std::floor(x/w+1e-14); double r=x-(double)q*w; while(r<0){ r+=w; q--; } while(r>=w){ r-=w; q++; } qval[i]=q; prefix_q[i+1]=prefix_q[i]+q; resid[i]=r; resid_sorted[i]=r; resid_order[i]=i; } std::sort(resid_sorted.begin(),resid_sorted.end()); std::sort(resid_order.begin(),resid_order.end(),[&](size_t a,size_t b){ if(resid[a]!=resid[b]) return resid[a]<resid[b]; return a<b; }); for(size_t i=0;i<n;i++) resid_rank[i]=(int)(std::lower_bound(resid_sorted.begin(),resid_sorted.end(),resid[i])-resid_sorted.begin()); }
    static ull u128_to_ull_checked(u128 x){ if(x>(u128)std::numeric_limits<ull>::max()) throw std::runtime_error("single layer count overflowed uint64"); return (ull)x; }
    u128 batch_sum_shift(const std::vector<double>& offsets,double T) const{ if(!base) return 0; const auto& S=base->sums; size_t mOff=(size_t)(std::upper_bound(offsets.begin(),offsets.end(),T+EPS)-offsets.begin()); Fenwick fw((int)resid_sorted.size()); size_t added=0; u128 ans=0; for(size_t rev=0;rev<mOff;++rev){ size_t qi=mOff-1-rev; double X=T-offsets[qi]; if(X<-EPS) continue; while(added<S.size() && S[added]<=X+EPS){ fw.add(resid_rank[added],1); added++; } if(added==0) continue; ull qx=(ull)std::floor(X/w+1e-14); double R=X-(double)qx*w; while(R<0){R+=w;qx--;} while(R>=w){R-=w;qx++;} size_t m=added; u128 term=(u128)m*(u128)(qx+1)-(u128)prefix_q[m]; int le_rank=(int)(std::upper_bound(resid_sorted.begin(),resid_sorted.end(),R+1e-14)-resid_sorted.begin())-1; int le=fw.sum_le_index(le_rank); ull gt=(ull)m-(ull)le; ans += term-(u128)gt; } return ans; }
    std::vector<ull> batch_values_shift(const std::vector<double>& offsets,double T,size_t mOff) const{ std::vector<ull> out(mOff,0); if(!base) return out; const auto& S=base->sums; Fenwick fw((int)resid_sorted.size()); size_t added=0; for(size_t rev=0;rev<mOff;++rev){ size_t qi=mOff-1-rev; double X=T-offsets[qi]; if(X<-EPS) continue; while(added<S.size() && S[added]<=X+EPS){ fw.add(resid_rank[added],1); added++; } if(added==0) continue; ull qx=(ull)std::floor(X/w+1e-14); double R=X-(double)qx*w; while(R<0){R+=w;qx--;} while(R>=w){R-=w;qx++;} size_t m=added; u128 term=(u128)m*(u128)(qx+1)-(u128)prefix_q[m]; int le_rank=(int)(std::upper_bound(resid_sorted.begin(),resid_sorted.end(),R+1e-14)-resid_sorted.begin())-1; int le=fw.sum_le_index(le_rank); ull gt=(ull)m-(ull)le; out[qi]=u128_to_ull_checked(term-(u128)gt); } return out; }
    void report_interval(double L,double H,std::vector<std::pair<size_t,int>>& out) const{ if(!base || H<-EPS) return; if(L<-EPS) L=-EPS; long long qlo=(long long)std::floor(L/w+1e-14), qhi=(long long)std::floor(H/w+1e-14); if(qhi<0) return; if(qlo<0) qlo=0; for(long long Q=qlo;Q<=qhi;++Q){ double low_r=(Q==qlo?L-(double)Q*w:-EPS); double high_r=(Q==qhi?H-(double)Q*w:w+EPS); if(high_r<-EPS||low_r>w+EPS) continue; if(low_r<-EPS) low_r=-EPS; if(high_r>w+EPS) high_r=w+EPS; auto it1=std::upper_bound(resid_order.begin(),resid_order.end(),low_r,[&](double val,size_t idx){return val<resid[idx];}); auto it2=std::upper_bound(resid_order.begin(),resid_order.end(),high_r+1e-14,[&](double val,size_t idx){return val<resid[idx];}); for(auto it=it1;it!=it2;++it){ size_t bi=*it; if(qval[bi]<=(ull)Q){ ull cu=(ull)Q-qval[bi]; if(cu<=(ull)std::numeric_limits<int>::max()){ double y=base->sums[bi]+(double)cu*w; if(y>L-EPS && y<=H+EPS) out.push_back({bi,(int)cu}); } } } } }
};
struct Candidate{ double sum; std::vector<int> exps; };

struct LayerRanker{
    std::vector<int> primes; std::vector<double> logs; double maxT; std::vector<int>Aidx,Baseidx; int outer=-1; Group A,Base; CompressedLayer layer; BuildStats astat,bst; double build_sec=0; double split_score=0;
    LayerRanker(std::vector<int> p,double maxT_):primes(std::move(p)),maxT(maxT_){ std::sort(primes.begin(),primes.end()); for(int q:primes) logs.push_back(std::log((double)q)); choose_split(); auto t0=std::chrono::high_resolution_clock::now(); A=gen_group_pointer_sorted(logs,Aidx,maxT,astat); Base=gen_group_pointer_sorted(logs,Baseidx,maxT,bst); layer.init(Base,logs[outer]); auto t1=std::chrono::high_resolution_clock::now(); build_sec=std::chrono::duration<double>(t1-t0).count(); }
    double approx_mask(int mask) const{ int k=(int)primes.size(),d=0; double prod=1.0; for(int i=0;i<k;i++) if(mask>>i&1){ d++; prod*=logs[i]; } if(d==0) return 1.0; double f=1.0; for(int j=2;j<=d;j++) f*=j; return std::pow(std::max(0.0,maxT),d)/(f*prod)+10.0*std::pow(std::max(1.0,maxT),std::max(0,d-1)); }
    void choose_split(){ int k=(int)primes.size(); if(k<3) throw std::runtime_error("layer ranker needs at least 3 primes"); int full=(1<<k)-1; double best=1e300; int bestA=0,bestOuter=-1; // enumerate A and one outer inside complement; complement-base is compressed side base
        for(int mask=1; mask<full; ++mask){ int ad=__builtin_popcount((unsigned)mask); int bdim=k-ad; if(bdim<2) continue; // need base + outer
            // keep A no larger than k-2, and not too tiny for high k unless it scores well.
            int bmask=full^mask; for(int o=0;o<k;o++) if(bmask>>o&1){ int basemask=bmask^(1<<o); double a=approx_mask(mask), base=approx_mask(basemask), bfull=approx_mask(bmask); // estimated cost: build A+base, ~14 count calls over A plus base additions. penalize giant base less than giant full B.
                    double score = a + base + 14.0*(0.20*a + 0.04*base) + 0.00000003*a*std::log(base+2.0); // heuristic
                    // prefer larger outer log because it compresses more layers
                    score *= (1.0 + 0.015*(logs[o]<*std::max_element(logs.begin(),logs.end())));
                    if(score<best){ best=score; bestA=mask; bestOuter=o; }
            }
        }
        split_score=best; int bmask=full^bestA; outer=bestOuter; for(int i=0;i<k;i++) if(bestA>>i&1) Aidx.push_back(i); for(int i=0;i<k;i++) if((bmask>>i&1) && i!=outer) Baseidx.push_back(i); if(Aidx.empty()||Baseidx.empty()||outer<0) throw std::runtime_error("bad split"); }
    u128 count_le(double T) const{ return layer.batch_sum_shift(A.sums,T); }
    std::vector<int> combine(Pack ap,Pack bp,int c) const{ std::vector<int> e(primes.size(),0); for(int i=0;i<(int)Aidx.size();i++) e[Aidx[i]]=pack_get(ap,i); for(int i=0;i<(int)Baseidx.size();i++) e[Baseidx[i]]=pack_get(bp,i); e[outer]=c; return e; }
    std::vector<Candidate> band_slow(double low,double high) const{ std::vector<Candidate> out; if(high<-EPS) return out; int cmax=(int)std::floor(high/logs[outer]+1e-14); for(int c=0;c<=cmax;c++){ double l=low-(double)c*logs[outer], u=high-(double)c*logs[outer]; if(u<-EPS) break; long long hi=(long long)Base.sums.size()-1, lo=(long long)Base.sums.size()-1; for(size_t ai=0; ai<A.sums.size(); ++ai){ double a=A.sums[ai]; if(a>u+EPS) break; double thu=u-a; while(hi>=0 && Base.sums[(size_t)hi]>thu+EPS) --hi; if(hi<0) break; double thl=l-a; while(lo>=0 && Base.sums[(size_t)lo]>thl+EPS) --lo; long long start=lo+1,end=hi; if(start<=end){ for(long long bj=start; bj<=end; ++bj){ double s=a+Base.sums[(size_t)bj]+(double)c*logs[outer]; if(s>low-EPS&&s<=high+EPS) out.push_back({s,combine(A.packs[ai],Base.packs[(size_t)bj],c)}); } } } } std::sort(out.begin(),out.end(),[](const Candidate&x,const Candidate&y){ if(x.sum!=y.sum) return x.sum<y.sum; return x.exps<y.exps; }); return out; }
    std::vector<Candidate> band(double low,double high) const{ std::vector<Candidate> out; if(high<-EPS) return out; size_t mA=(size_t)(std::upper_bound(A.sums.begin(),A.sums.end(),high+EPS)-A.sums.begin()); auto ch=layer.batch_values_shift(A.sums,high,mA); auto cl=layer.batch_values_shift(A.sums,low,mA); ull expected=0; std::vector<size_t> active; active.reserve(1024); for(size_t i=0;i<mA;i++){ if(ch[i]>cl[i]){ expected += ch[i]-cl[i]; active.push_back(i); } } out.reserve((size_t)std::min<ull>(expected,10000000ULL)); std::vector<std::pair<size_t,int>> reps; reps.reserve(128); for(size_t ai: active){ reps.clear(); double xl=low-A.sums[ai], xh=high-A.sums[ai]; layer.report_interval(xl,xh,reps); for(auto &bc:reps){ size_t bj=bc.first; int c=bc.second; double s=A.sums[ai]+Base.sums[bj]+(double)c*logs[outer]; if(s>low-EPS&&s<=high+EPS) out.push_back({s,combine(A.packs[ai],Base.packs[bj],c)}); } } if(out.size()!=(size_t)expected) return band_slow(low,high); std::sort(out.begin(),out.end(),[](const Candidate&x,const Candidate&y){ if(x.sum!=y.sum) return x.sum<y.sum; return x.exps<y.exps; }); return out; }
};
struct RankResult{ double seconds=0,build=0,count_phase=0,band_phase=0,exact=0,leading_seed=0,analytic_seed=0; bool analytic_bracket=false; std::vector<int> exps; int digits=0; size_t A=0,Base=0,band=0; int calls=0; u128 rank_gap=0; std::vector<int>Aidx,Baseidx; int outer=-1; };
static RankResult nth_layer(std::vector<int> primes,ull n,ull target_gap=50000ULL){ auto t0=std::chrono::high_resolution_clock::now(); double lead=leading_est(primes,n); double aest=asymptotic_est(primes,n); bool use_analytic=(n>=1000000ULL && aest>0.0 && std::isfinite(aest) && aest>0.5*lead && aest<1.5*lead); ull effective_gap=use_analytic?std::min<ull>(target_gap,5000ULL):target_gap; double est=use_analytic?aest:lead; double maxT=std::max(1e-9,std::max(lead*1.02,est*1.01)+1e-8); std::unique_ptr<LayerRanker> R(new LayerRanker(primes,maxT)); int calls=0; auto count=[&](double T){ calls++; return R->count_le(T); }; double lo=use_analytic?std::max(0.0,est*0.985):std::max(0.0,lead*0.70), hi=R->maxT; if(lo>=hi) lo=std::max(0.0,hi*0.75); u128 c_hi=count(hi); while(c_hi<(u128)n){ maxT*=1.25; R.reset(new LayerRanker(primes,maxT)); hi=maxT; c_hi=count(hi); } u128 c_lo=count(lo); while(c_lo>=(u128)n){ hi=lo; c_hi=c_lo; lo*=0.75; c_lo=count(lo); }
    auto tc0=std::chrono::high_resolution_clock::now(); int stagnant=0; for(int iter=0;iter<44;iter++){ u128 gap=c_hi-c_lo; if(gap<=(u128)effective_gap) break; double t; long double denom=(long double)(c_hi-c_lo); if(denom<=0) t=(lo+hi)*0.5; else{ long double frac=((long double)n-(long double)c_lo)/denom; if(!(frac>0.0L&&frac<1.0L)) frac=0.5L; if(frac<0.01L) frac=0.01L; if(frac>0.99L) frac=0.99L; t=lo+(double)frac*(hi-lo); } if(t<=lo+1e-14*(1+fabs(lo))||t>=hi-1e-14*(1+fabs(hi))||stagnant>=4){ t=(lo+hi)*0.5; stagnant=0; } u128 c=count(t); u128 oldgap=gap; if(c>=(u128)n){ hi=t; c_hi=c; } else { lo=t; c_lo=c; } u128 newgap=c_hi-c_lo; if(newgap>oldgap*95/100) stagnant++; else stagnant=0; }
    for(int iter=0;iter<50 && c_hi-c_lo>(u128)effective_gap;iter++){ double mid=(lo+hi)*0.5; u128 c=count(mid); if(c>=(u128)n){ hi=mid; c_hi=c; } else { lo=mid; c_lo=c; } }
    auto tc1=std::chrono::high_resolution_clock::now(); auto tb0=std::chrono::high_resolution_clock::now(); std::vector<Candidate> cands; u128 below=0; long long off=-1; double L=lo,H=hi; for(int attempt=0;attempt<8;attempt++){ below=count(L); cands=R->band(L,H); if(below<=(u128)n-1){ u128 off128=(u128)n-below-1; if(off128<=(u128)std::numeric_limits<long long>::max()) off=(long long)off128; } if(off>=0 && (size_t)off<cands.size()) break; double w=H-L; if(w<=0) w=1e-10*(1+fabs(H)); L=std::max(-1e-9,L-w); H=H+w; }
    auto tb1=std::chrono::high_resolution_clock::now(); if(!(off>=0 && (size_t)off<cands.size())) throw std::runtime_error("band failure below="+u128str(below)+" cands="+std::to_string(cands.size())+" off="+std::to_string(off)); auto te0=std::chrono::high_resolution_clock::now(); struct EC{Candidate c; cpp_int val;}; std::vector<EC> exact; exact.reserve(cands.size()); for(auto &c:cands) exact.push_back({c,value_from_exps(primes,c.exps)}); std::sort(exact.begin(),exact.end(),[](const EC&a,const EC&b){return a.val<b.val;}); auto chosen=exact[(size_t)off]; auto te1=std::chrono::high_resolution_clock::now(); auto t1=std::chrono::high_resolution_clock::now(); RankResult r; r.seconds=std::chrono::duration<double>(t1-t0).count(); r.build=R->build_sec; r.count_phase=std::chrono::duration<double>(tc1-tc0).count(); r.band_phase=std::chrono::duration<double>(tb1-tb0).count(); r.exact=std::chrono::duration<double>(te1-te0).count(); r.leading_seed=lead; r.analytic_seed=aest; r.analytic_bracket=use_analytic; r.exps=chosen.c.exps; r.digits=digits10(chosen.val); r.A=R->A.sums.size(); r.Base=R->Base.sums.size(); r.band=cands.size(); r.calls=calls; r.rank_gap=c_hi-c_lo; r.Aidx=R->Aidx; r.Baseidx=R->Baseidx; r.outer=R->outer; return r; }
static std::vector<int> first_primes(int k){ int base[]={2,3,5,7,11,13,17,19,23,29,31,37,41,43,47}; return std::vector<int>(base,base+k); }
int main(int argc,char**argv){ std::cout.setf(std::ios::fixed); std::cout<<std::setprecision(6); try{ if(argc>=2 && std::string(argv[1])=="nth"){ auto P=parse_csv(argv[2]); ull n=std::stoull(argv[3]); ull gap=argc>=5?std::stoull(argv[4]):50000ULL; auto r=nth_layer(P,n,gap); std::cout<<"layer_compressed P="<<join_primes(P)<<" k="<<P.size()<<" N="<<n<<" seconds="<<r.seconds<<" build="<<r.build<<" count_phase="<<r.count_phase<<" band_phase="<<r.band_phase<<" exact="<<r.exact<<" calls="<<r.calls<<" rank_gap="<<u128str(r.rank_gap)<<" leading_seed="<<r.leading_seed<<" analytic_seed="<<r.analytic_seed<<" analytic_bracket="<<(r.analytic_bracket?"true":"false")<<" A="<<r.A<<" Base="<<r.Base<<" band="<<r.band<<" exps="<<join_vec(r.exps)<<" digits="<<r.digits<<" splitA="<<idx_to_primes(r.Aidx,P)<<" splitBase="<<idx_to_primes(r.Baseidx,P)<<" outer="<<P[r.outer]<<"\n"; return 0; } for(int k:{5,6,8}){ auto P=first_primes(k); auto r=nth_layer(P,1000000000000ULL,50000ULL); std::cout<<"layer_compressed P="<<join_primes(P)<<" k="<<P.size()<<" N=1000000000000 seconds="<<r.seconds<<" build="<<r.build<<" count_phase="<<r.count_phase<<" band_phase="<<r.band_phase<<" exact="<<r.exact<<" calls="<<r.calls<<" rank_gap="<<u128str(r.rank_gap)<<" leading_seed="<<r.leading_seed<<" analytic_seed="<<r.analytic_seed<<" analytic_bracket="<<(r.analytic_bracket?"true":"false")<<" A="<<r.A<<" Base="<<r.Base<<" band="<<r.band<<" exps="<<join_vec(r.exps)<<" digits="<<r.digits<<" splitA="<<idx_to_primes(r.Aidx,P)<<" splitBase="<<idx_to_primes(r.Baseidx,P)<<" outer="<<P[r.outer]<<"\n"; } }catch(const std::exception&e){ std::cerr<<"error: "<<e.what()<<"\n"; return 1; } }
