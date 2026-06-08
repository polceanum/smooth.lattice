#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <iomanip>
#include <array>
#include <limits>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>
#include <boost/multiprecision/cpp_int.hpp>

using ull = unsigned long long;
using u128 = __uint128_t;
using boost::multiprecision::cpp_int;
static constexpr long double EPS = 1e-18L;

static std::vector<int> parse_csv(const std::string& s){
    std::vector<int> out; std::stringstream ss(s); std::string item;
    while(std::getline(ss,item,',')) if(!item.empty()) out.push_back(std::stoi(item));
    std::sort(out.begin(), out.end()); return out;
}
static std::string join_vec(const std::vector<int>& v){ std::ostringstream os; os<<"["; for(size_t i=0;i<v.size();++i){ if(i) os<<","; os<<v[i]; } os<<"]"; return os.str(); }
static std::string join_primes(const std::vector<int>& v){ std::ostringstream os; os<<"("; for(size_t i=0;i<v.size();++i){ if(i) os<<","; os<<v[i]; } os<<")"; return os.str(); }
static cpp_int ipow(int base, int exp){ cpp_int r=1,b=base; while(exp){ if(exp&1) r*=b; exp>>=1; if(exp) b*=b; } return r; }
static cpp_int value_from_exps(const std::vector<int>& primes, const std::vector<int>& exps){ cpp_int r=1; for(size_t i=0;i<primes.size();++i) if(exps[i]) r*=ipow(primes[i], exps[i]); return r; }
static int digits10(const cpp_int& x){ std::stringstream ss; ss<<x; return (int)ss.str().size(); }
static unsigned long long factorial(int k){ unsigned long long r=1; for(int i=2;i<=k;i++) r*=i; return r; }
static long double leading_est(const std::vector<int>& primes, ull n){ long double prod=1.0L; for(int p:primes) prod*=logl((long double)p); int k=(int)primes.size(); return powl((long double)n*(long double)factorial(k)*prod, 1.0L/(long double)k); }
static std::string u128str(u128 x){ if(x==0) return "0"; std::string s; while(x){ int d=x%10; s.push_back('0'+d); x/=10; } std::reverse(s.begin(),s.end()); return s; }

// Beatty/floor-sum core: sum_{i=0}^{n-1} floor(a*i+b). Depth is O(log n) from continued fractions.
static u128 floor_sum_rec(ull n, long double a, long double b){
    if(n == 0) return 0;
    if(a < 1e-30L) return (u128)n * (ull)floorl(b + EPS);
    u128 ans = 0;
    long double fald = floorl(a + EPS);
    if(fald >= 1.0L){ ull fa=(ull)fald; ans += (u128)fa*(u128)n*(u128)(n-1)/2; a -= (long double)fa; }
    long double fbld = floorl(b + EPS);
    if(fbld >= 1.0L){ ull fb=(ull)fbld; ans += (u128)fb*(u128)n; b -= (long double)fb; }
    if(a < 1e-30L) return ans;
    ull m = (ull)floorl(a*(long double)n + b + EPS);
    if(m == 0) return ans;
    ans += (u128)(n-1) * (u128)m;
    ans -= floor_sum_rec(m, 1.0L/a, (1.0L-b)/a);
    return ans;
}

static u128 count2_loop(long double T, long double w0, long double w1){
    if(T < -1e-18L) return 0;
    // Loop over larger weight to keep loops smaller.
    long double a = std::min(w0,w1), b = std::max(w0,w1);
    ull m = (ull)floorl(T/b + EPS);
    u128 total = 0;
    for(ull j=0;j<=m;j++){
        long double rem = T - (long double)j*b;
        if(rem >= -EPS) total += (ull)floorl(rem/a + EPS) + 1ULL;
    }
    return total;
}
static u128 count2_beatty(long double T, long double w0, long double w1){
    if(T < -1e-18L) return 0;
    long double a = std::min(w0,w1), b = std::max(w0,w1);
    ull m = (ull)floorl(T/b + EPS);
    long double rem = T - (long double)m*b; // in [0,b)
    // Sum_{j=0}^m ( floor((rem + j*b)/a) + 1 )
    return (u128)(m+1) + floor_sum_rec(m+1, b/a, rem/a);
}

struct Candidate{ long double logsum; std::vector<int> exps; cpp_int val; };

struct ThreeBeattyRanker{
    std::vector<int> primes;
    std::vector<long double> logs;
    int outer; std::vector<int> inner;
    ThreeBeattyRanker(std::vector<int> p): primes(std::move(p)){
        std::sort(primes.begin(), primes.end());
        if(primes.size()!=3) throw std::runtime_error("ThreeBeattyRanker requires exactly 3 primes");
        for(int q:primes) logs.push_back(logl((long double)q));
        outer = 0;
        for(int i=1;i<3;i++) if(logs[i] > logs[outer]) outer = i;
        for(int i=0;i<3;i++) if(i!=outer) inner.push_back(i);
    }
    u128 count_le(long double T) const{
        if(T < -1e-18L) return 0;
        long double wo = logs[outer];
        ull mo = (ull)floorl(T/wo + EPS);
        u128 total=0;
        for(ull e=0;e<=mo;e++){
            total += count2_beatty(T - (long double)e*wo, logs[inner[0]], logs[inner[1]]);
        }
        return total;
    }
    void enum2_strip(long double L, long double H, int eouter, std::vector<Candidate>& out) const{
        int i0=inner[0], i1=inner[1];
        long double w0=logs[i0], w1=logs[i1];
        // Loop over larger inner weight to reduce iterations.
        int loopIdx=i0, solveIdx=i1;
        long double wLoop=w0, wSolve=w1;
        if(w0 < w1){ loopIdx=i1; solveIdx=i0; wLoop=w1; wSolve=w0; }
        if(H < -EPS) return;
        ull m=(ull)floorl(H/wLoop + EPS);
        for(ull eloop=0; eloop<=m; ++eloop){
            long double base = (long double)eloop*wLoop;
            long double lowNeed = (L - base)/wSolve;
            long double highNeed = (H - base)/wSolve;
            long long amin = (long long)floorl(lowNeed + EPS) + 1; // strict > L
            long long amax = (long long)floorl(highNeed + EPS);     // <= H
            if(amin < 0) amin = 0;
            if(amax < amin) continue;
            for(long long esolve=amin; esolve<=amax; ++esolve){
                std::vector<int> e(3,0);
                e[outer]=eouter; e[loopIdx]=(int)eloop; e[solveIdx]=(int)esolve;
                long double pair_s = logs[loopIdx]*(long double)eloop + logs[solveIdx]*(long double)esolve;
                long double s = pair_s + logs[outer]*(long double)eouter;
                if(pair_s > L - 1e-15L && pair_s <= H + 1e-15L){
                    cpp_int v=value_from_exps(primes,e);
                    out.push_back({s,e,v});
                }
            }
        }
    }
    std::vector<Candidate> band(long double L, long double H) const{
        std::vector<Candidate> out;
        if(H < -EPS) return out;
        long double wo=logs[outer];
        ull mo=(ull)floorl(H/wo + EPS);
        for(ull e=0;e<=mo;e++){
            long double L2=L-(long double)e*wo, H2=H-(long double)e*wo;
            u128 c = count2_beatty(H2, logs[inner[0]], logs[inner[1]]) - count2_beatty(L2, logs[inner[0]], logs[inner[1]]);
            if(c>0) enum2_strip(L2,H2,(int)e,out);
        }
        std::sort(out.begin(), out.end(), [](const Candidate& a, const Candidate& b){
            if(a.val != b.val) return a.val < b.val;
            return a.exps < b.exps;
        });
        return out;
    }
};

struct Result{ double seconds; std::vector<int> exps; int digits; size_t band; int calls; long double logsum; };
static Result nth3_beatty(std::vector<int> primes, ull n){
    auto t0=std::chrono::high_resolution_clock::now();
    ThreeBeattyRanker R(primes);
    long double est=leading_est(primes,n);
    long double lo=std::max((long double)0.0, est*0.80L), hi=std::max((long double)1e-9, est*1.20L + 1e-8L);
    int calls=0;
    auto count=[&](long double T){ calls++; return R.count_le(T); };
    while(count(hi) < (u128)n) hi*=1.25L;
    while(count(lo) >= (u128)n){ hi=lo; lo*=0.75L; }
    for(int i=0;i<70;i++){
        long double mid=(lo+hi)/2;
        if(count(mid) >= (u128)n) hi=mid; else lo=mid;
    }
    long double center=hi;
    long double width=1e-8L;
    std::vector<Candidate> cands;
    u128 below=0; long long off=-1;
    for(int attempt=0; attempt<16; ++attempt){
        long double L=center-width, H=center+width;
        below=count(L);
        cands=R.band(L,H);
        if((u128)n > below){
            u128 off128 = (u128)n - below - 1;
            if(off128 <= (u128)std::numeric_limits<long long>::max()) off=(long long)off128;
            else off=-1;
        } else off=-1;
        if(off>=0 && (size_t)off<cands.size()) break;
        width *= 5.0L;
    }
    if(!(off>=0 && (size_t)off<cands.size())){
        std::ostringstream os; os << "band failure: cands="<<cands.size()<<" off="<<off<<" below="<<u128str(below);
        throw std::runtime_error(os.str());
    }
    Candidate chosen=cands[(size_t)off];
    auto t1=std::chrono::high_resolution_clock::now();
    return {std::chrono::duration<double>(t1-t0).count(), chosen.exps, digits10(chosen.val), cands.size(), calls, chosen.logsum};
}

// Simple pointer/log baseline with exponent vectors, for validation at moderate N.
static std::vector<int> pointer_exps_nth(const std::vector<int>& primes, ull n){
    int k=(int)primes.size(); std::vector<long double> L(k); for(int i=0;i<k;i++) L[i]=logl((long double)primes[i]);
    std::vector<long double> H((size_t)n); H[0]=0.0L;
    std::vector<std::array<int,3>> E((size_t)n); E[0]={0,0,0};
    std::vector<ull> idx(k,0);
    for(ull t=1;t<n;t++){
        long double nxt=H[(size_t)idx[0]]+L[0]; int src=0;
        for(int j=1;j<k;j++){ long double v=H[(size_t)idx[j]]+L[j]; if(v<nxt){ nxt=v; src=j; } }
        H[(size_t)t]=nxt; E[(size_t)t]=E[(size_t)idx[src]]; E[(size_t)t][src]++;
        for(int j=0;j<k;j++) while(idx[j]<t && H[(size_t)idx[j]]+L[j] <= nxt + 1e-15L) idx[j]++;
    }
    return {E[(size_t)n-1][0],E[(size_t)n-1][1],E[(size_t)n-1][2]};
}

int main(int argc, char** argv){
    std::cout.setf(std::ios::fixed); std::cout<<std::setprecision(6);
    if(argc>=3 && std::string(argv[1])=="countcheck"){
        auto P=parse_csv(argv[2]); long double T=strtold(argv[3],nullptr); ThreeBeattyRanker R(P);
        std::cout << "count3 P="<<join_primes(P)<<" T="<<(double)T<<" count="<<u128str(R.count_le(T))<<"\n";
        return 0;
    }
    if(argc>=4 && std::string(argv[1])=="nth"){
        auto P=parse_csv(argv[2]); ull n=std::stoull(argv[3]);
        auto r=nth3_beatty(P,n);
        std::cout << "beatty3_rank P="<<join_primes(P)<<" N="<<n<<" seconds="<<r.seconds<<" exps="<<join_vec(r.exps)<<" digits="<<r.digits<<" band="<<r.band<<" calls="<<r.calls<<" log="<<std::setprecision(12)<<(double)r.logsum<<std::setprecision(6)<<"\n";
        return 0;
    }
    if(argc>=4 && std::string(argv[1])=="validate"){
        auto P=parse_csv(argv[2]); ull n=std::stoull(argv[3]);
        auto r=nth3_beatty(P,n); auto e=pointer_exps_nth(P,n);
        std::cout << "validate P="<<join_primes(P)<<" N="<<n<<" beatty="<<join_vec(r.exps)<<" pointer="<<join_vec(e)<<" match="<<(r.exps==e?"true":"false")<<" seconds="<<r.seconds<<"\n";
        return 0;
    }
    // self-test count2 against loops for a few cases
    std::vector<std::pair<int,int>> pairs={{2,3},{2,5},{3,5},{5,7}};
    for(auto [p,q]: pairs){
        long double a=logl((long double)p), b=logl((long double)q);
        for(long double T: {1.0L,5.0L,20.0L,100.0L,1000.0L}){
            auto x=count2_beatty(T,a,b), y=count2_loop(T,a,b);
            std::cout << "count2_check ("<<p<<","<<q<<") T="<<(double)T<<" beatty="<<u128str(x)<<" loop="<<u128str(y)<<" match="<<(x==y?"true":"false")<<"\n";
        }
    }
    std::vector<std::vector<int>> Ps={{2,3,5},{2,3,7},{2,5,7},{3,5,7}};
    std::vector<ull> Ns={1000000ULL,1000000000ULL,1000000000000ULL};
    for(auto P:Ps){
        for(ull n:Ns){
            auto r=nth3_beatty(P,n);
            std::cout << "beatty3_rank P="<<join_primes(P)<<" N="<<n<<" seconds="<<r.seconds<<" exps="<<join_vec(r.exps)<<" digits="<<r.digits<<" band="<<r.band<<" calls="<<r.calls<<"\n";
        }
    }
    for(auto P:Ps){ for(ull n: {1000ULL,1000000ULL}){
        auto r=nth3_beatty(P,n); auto e=pointer_exps_nth(P,n);
        std::cout << "validate P="<<join_primes(P)<<" N="<<n<<" match="<<(r.exps==e?"true":"false")<<" beatty="<<join_vec(r.exps)<<" pointer="<<join_vec(e)<<"\n";
    }}
    return 0;
}
