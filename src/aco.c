#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "aco.h"
#include "eval.h"
#include "util.h"

#include <stdio.h>
#define TRACE(fmt, ...) do { fprintf(stderr, "[ACO] " fmt "\n", ##__VA_ARGS__); fflush(stderr); } while(0)
#define SOFT_ENABLE_CC 1
#define SOFT_PASSES 3


// -------------------- helpers & utilities --------------------

static inline int day_of(const Instance* I, int t){ return t / I->periods_per_day; }
static inline int per_of(const Instance* I, int t){ return t % I->periods_per_day; }

static int L2C(const Instance* I, int L){
    return I->lecture_to_course[L];
}

static double** new_matrix(int R, int C, double v){
    double** M = (double**)malloc(R * sizeof(*M));
    for(int r=0;r<R;++r){
        M[r] = (double*)malloc(C*sizeof(**M));
        for(int c=0;c<C;++c) M[r][c]=v;
    }
    return M;
}
static void free_matrix(double** M, int R){
    for(int r=0;r<R;++r) free(M[r]);
    free(M);
}

// Heuristic: prefer times where a fitting room exists (tight fit),
// otherwise penalize by overflow; add small pressure for conflict degree.
// Heuristic: capacity-first (tight fit > loose fit > overflow), plus small conflict pressure.
static double heuristic_eta(Instance* I, int L, int t){
    int c = L2C(I, L);
    if(!I->feasible[c][t]) return 1e-12;

    int need = I->courses[c].students;

    // scan rooms at this timeslot
    int fits = 0;
    int slack_best = 1000000000;   // smallest non-negative slack (cap - need)
    int min_over   = 1000000000;   // smallest overflow (need - cap), positive only

    for(int r=0; r<I->num_rooms; ++r){
        int cap   = I->rooms[r].capacity;
        int slack = cap - need;           // >=0 means it fits with 'slack' seats
        int over  = need - cap;           // >0 means overflow by 'over' seats

        if(slack >= 0){                   // fits
            fits = 1;
            if(slack < slack_best) slack_best = slack;
        } else {                          // overflow
            if(over < min_over) min_over = over;
        }
    }

    // capacity-driven base score
    double base;
    if(fits){
        // Prefer tighter fit (slack close to 0 -> higher score)
        // Tunables: 6.0 (bias for any fit), 3.0 (how much we reward tighter fits)
        base = 6.0 + 3.0 / (1.0 + (double)slack_best);
    }else{
        // No room fits: penalize proportional to overflow (scale 0.02 is tunable)
        double of = (double)(min_over > 0 ? min_over : 0);
        base = 1.0 / (1.0 + 0.02 * of);
    }
    if(base < 1e-9) base = 1e-9;

    // small static pressure for conflict degree (more conflicts -> slightly higher score)
    int deg = 0;
    for(int cc=0; cc<I->num_courses; ++cc) if(I->conflict[c][cc]) ++deg;

    return base * (1.0 + 0.01 * deg);
}

// Best-fit room that avoids room-time collisions; prefer smallest that fits.
static int choose_room(const Instance* I, int c, int t, int occ_room_time[][MAX_TIMESLOTS]){
    int best = -1, best_cap = 1000000000, best_over = 1000000000;
    for(int r=0; r<I->num_rooms; ++r){
        if(occ_room_time[r][t] > 0) continue; // already used at (r,t)
        int cap  = I->rooms[r].capacity;
        int need = I->courses[c].students;
        int over = need - cap; // <= 0 means fits
        if(over <= 0){
            if(cap < best_cap){ best = r; best_cap = cap; best_over = over; }
        }else if(best == -1 && over < best_over){
            best = r; best_over = over;
        }
    }
    if(best == -1) best = 0; // fallback
    return best;
}

// Top-p% selective timeslots by heuristic (among feasible)
static int build_selective_candidates(Instance* I, int L, double selective_p, int* cand, int max){
    double score[MAX_TIMESLOTS];
    int    list [MAX_TIMESLOTS];
    int n=0;

    for(int t=0; t<I->num_timeslots; ++t){
        if(I->feasible[L2C(I,L)][t]){
            list[n]=t; score[n]=heuristic_eta(I,L,t); ++n;
        }
    }
    if(n==0) return 0;

    int k = (int)ceil(selective_p * n);
    if(k<1) k=1; if(k>max) k=max;

    // partial selection for top-k
    for(int i=0;i<k;++i){
        int best=i;
        for(int j=i+1;j<n;++j) if(score[j] > score[best]) best=j;
        int tt=list[i]; list[i]=list[best]; list[best]=tt;
        double ss=score[i]; score[i]=score[best]; score[best]=ss;
        cand[i]=list[i];
    }
    return k;
}

// Difficulty = lectures with few feasible slots, many students, high conflict degree
static void compute_difficulty(const Instance* I, double* diff){ // size = total_lectures
    for(int L=0; L<I->total_lectures; ++L){
        int c = I->lecture_to_course[L];
        // feasible slot count
        int feas=0; for(int t=0;t<I->num_timeslots;++t) if(I->feasible[c][t]) ++feas;
        if(feas==0) feas=1; // avoid div-by-zero, treat as hardest

        // conflict degree (course-level)
        int deg=0; for(int cc=0; cc<I->num_courses; ++cc) if(I->conflict[c][cc]) ++deg;

        // students
        int stud = I->courses[c].students;

        // smaller 'feas' -> bigger difficulty, so weight 1/feas
        diff[L] = 1000.0*(1.0/feas) + 0.1*deg + 0.01*stud;
    }
}

// Make an index order of lectures sorted by descending difficulty
static void make_order_desc(const Instance* I, int* order, double* diff){ // arrays len = total_lectures
    for(int L=0; L<I->total_lectures; ++L) order[L]=L;
    // simple selection sort (OK for these sizes)
    for(int i=0;i<I->total_lectures;++i){
        int best=i;
        for(int j=i+1;j<I->total_lectures;++j)
            if(diff[order[j]] > diff[order[best]]) best=j;
        int tmp = order[i]; order[i]=order[best]; order[best]=tmp;
    }
}

// -------------------- repairs --------------------

// Greedy multi-pass repair:
// (a) move lectures from unavailable slots
// (b) deconflict room-time collisions
// (c) spread days to meet MinWorkingDays
// (d) fix teacher/curriculum clashes per timeslot
// (e) place any unassigned lectures (-1) to a feasible, non-clashing slot/room
static void repair_basic(Instance* I, Timetable* S){
    // run a few passes so changes can cascade
    for(int pass=0; pass<3; ++pass){
        // a) fix infeasible placements
        for(int L=0; L<I->total_lectures; ++L){
            int t = S->lec_to_timeslot[L];
            if(t < 0) continue;
            int c = I->lecture_to_course[L];
            if(!I->feasible[c][t]){
                for(int tt=0; tt<I->num_timeslots; ++tt){
                    if(I->feasible[c][tt]){ S->lec_to_timeslot[L]=tt; break; }
                }
            }
        }

        // b) resolve room-time collisions
        static int occ[MAX_ROOMS][MAX_TIMESLOTS];
        memset(occ,0,sizeof(occ));
        for(int L=0; L<I->total_lectures; ++L){
            int t=S->lec_to_timeslot[L]; if(t<0) continue;
            int r=S->lec_to_room[L];
            if(r>=0 && r<I->num_rooms) occ[r][t]++;
        }
        for(int L=0; L<I->total_lectures; ++L){
            int t=S->lec_to_timeslot[L]; if(t<0) continue;
            int r=S->lec_to_room[L];
            if(r<0 || r>=I->num_rooms) continue;
            if(occ[r][t] > 1){
                int c = I->lecture_to_course[L];
                // try same t, different room
                int moved=0;
                for(int rr=0; rr<I->num_rooms; ++rr){
                    if(rr==r) continue;
                    if(occ[rr][t]==0){
                        S->lec_to_room[L]=rr; occ[rr][t]++; occ[r][t]--; moved=1; break;
                    }
                }
                if(!moved){
                    // try feasible alt t with empty room (least overflow)
                    int best_rr=-1, best_tt=-1, best_over=1000000000;
                    for(int tt=0; tt<I->num_timeslots; ++tt){
                        if(!I->feasible[c][tt]) continue;
                        // avoid clashes at tt
                        int clash=0;
                        for(int L2=0; L2<I->total_lectures; ++L2)
                            if(S->lec_to_timeslot[L2]==tt && I->conflict[c][ I->lecture_to_course[L2] ]){ clash=1; break; }
                        if(clash) continue;
                        for(int rr=0; rr<I->num_rooms; ++rr){
                            if(occ[rr][tt]>0) continue;
                            int over = I->courses[c].students - I->rooms[rr].capacity;
                            if(over < best_over){ best_over = over; best_rr = rr; best_tt = tt; }
                        }
                    }
                    if(best_rr>=0){
                        occ[r][t]--;
                        S->lec_to_timeslot[L]=best_tt;
                        S->lec_to_room[L]=best_rr;
                        occ[best_rr][best_tt]++;
                    }
                }
            }
        }

        // c) improve MinWorkingDays by spreading to new days when short
        for (int c = 0; c < I->num_courses; ++c) {
            int start = I->course_first_lecture[c];
            int cnt   = I->course_lecture_count[c];
            int day_mask = 0;
            for (int k = 0; k < cnt; ++k) {
                int L = start + k, t = S->lec_to_timeslot[L];
                if (t >= 0) day_mask |= (1 << (t / I->periods_per_day));
            }
            int used = 0; for (int d = 0; d < I->days; ++d) if (day_mask & (1 << d)) ++used;
            int need = I->courses[c].min_working_days;
            if (used >= need) continue;

            // move one lecture to a new day
            for (int k = 0; k < cnt && used < need; ++k) {
                int L = start + k;
                int t0 = S->lec_to_timeslot[L];
                if (t0 < 0) continue;
                for (int d = 0; d < I->days; ++d) {
                    if (day_mask & (1 << d)) continue; // already used
                    for (int p = 0; p < I->periods_per_day; ++p) {
                        int tt = d * I->periods_per_day + p;
                        if (!I->feasible[c][tt]) continue;

                        // ensure no clash at tt
                        int clash=0;
                        for(int L2=0; L2<I->total_lectures; ++L2)
                            if(S->lec_to_timeslot[L2]==tt && I->conflict[c][ I->lecture_to_course[L2] ]){ clash=1; break; }
                        if(clash) continue;

                        // choose any empty room (prefer least overflow)
                        int best_rr = -1, best_over = 1000000000;
                        for (int rr = 0; rr < I->num_rooms; ++rr) {
                            // room free at tt?
                            int taken=0;
                            for(int L2=0; L2<I->total_lectures; ++L2)
                                if(S->lec_to_timeslot[L2]==tt && S->lec_to_room[L2]==rr){ taken=1; break; }
                            if(taken) continue;
                            int over = I->courses[c].students - I->rooms[rr].capacity;
                            if (over < best_over) { best_over = over; best_rr = rr; }
                        }
                        if (best_rr >= 0) {
                            S->lec_to_timeslot[L] = tt;
                            S->lec_to_room[L]     = best_rr;
                            day_mask |= (1 << d);
                            ++used;
                            goto next_course_mwd;
                        }
                    }
                }
            }
            next_course_mwd: ;
        }

        // d) fix teacher/curriculum clashes at each timeslot (greedy local move)
        for(int t=0; t<I->num_timeslots; ++t){
            int *Ls = (int*)malloc(sizeof(int)*I->total_lectures);
            int nL=0;
            for(int L=0; L<I->total_lectures; ++L)
                if(S->lec_to_timeslot[L]==t) Ls[nL++]=L;

            for(int i=0;i<nL;i++){
                int Li = Ls[i], ci = I->lecture_to_course[Li];
                for(int j=i+1;j<nL;j++){
                    int Lj = Ls[j], cj = I->lecture_to_course[Lj];
                    if(!I->conflict[ci][cj]) continue;

                    int moveL = (I->courses[ci].students <= I->courses[cj].students) ? Li : Lj;
                    int cm    = I->lecture_to_course[moveL];

                    int best_tt=-1, best_rr=-1, best_over=1000000000;
                    for(int dt=-2; dt<=2; ++dt){
                        int tt = t + dt;
                        if(tt<0 || tt>=I->num_timeslots) continue;
                        if(!I->feasible[cm][tt]) continue;

                        int clash2=0;
                        for(int L2=0; L2<I->total_lectures; ++L2)
                            if(S->lec_to_timeslot[L2]==tt && I->conflict[cm][ I->lecture_to_course[L2] ]){ clash2=1; break; }
                        if(clash2) continue;

                        for(int rr=0; rr<I->num_rooms; ++rr){
                            int taken=0;
                            for(int L2=0; L2<I->total_lectures; ++L2)
                                if(S->lec_to_timeslot[L2]==tt && S->lec_to_room[L2]==rr){ taken=1; break; }
                            if(taken) continue;

                            int over = I->courses[cm].students - I->rooms[rr].capacity;
                            if(over < best_over){ best_over=over; best_rr=rr; best_tt=tt; }
                        }
                    }
                    if(best_tt>=0){
                        S->lec_to_timeslot[moveL]=best_tt;
                        S->lec_to_room[moveL]=best_rr;
                    }
                }
            }
            free(Ls);
        }

        // e) place any unassigned lectures
        for(int L=0; L<I->total_lectures; ++L){
            if(S->lec_to_timeslot[L] >= 0) continue;
            int c = I->lecture_to_course[L];
            int best_tt=-1, best_rr=-1, best_over=1000000000;
            for(int tt=0; tt<I->num_timeslots; ++tt){
                if(!I->feasible[c][tt]) continue;
                // avoid clashes at tt
                int clash=0;
                for(int L2=0; L2<I->total_lectures; ++L2)
                    if(S->lec_to_timeslot[L2]==tt && I->conflict[c][ I->lecture_to_course[L2] ]){ clash=1; break; }
                if(clash) continue;
                // any free room?
                for(int rr=0; rr<I->num_rooms; ++rr){
                    int taken=0;
                    for(int L2=0; L2<I->total_lectures; ++L2)
                        if(S->lec_to_timeslot[L2]==tt && S->lec_to_room[L2]==rr){ taken=1; break; }
                    if(taken) continue;
                    int over = I->courses[c].students - I->rooms[rr].capacity;
                    if(over < best_over){ best_over=over; best_rr=rr; best_tt=tt; }
                }
            }
            if(best_tt>=0){
                S->lec_to_timeslot[L]=best_tt;
                S->lec_to_room[L]=best_rr;
            }
        }
    } // end passes
}

// -------------------- candidate weighting with current occupancy --------------------

// candidate weight with current occupancy/clash map and capacity awareness
// candidate weight with current occupancy/clash map and capacity awareness
static double weight_for(ACOParams P, Instance* I, int L, int t,
                         double tauLt,
                         int (*occ_room_time)[MAX_TIMESLOTS],
                         int *slot_used)
{
    int c = I->lecture_to_course[L];
    if (t < 0 || t >= I->num_timeslots) return 0.0;
    if (!I->feasible[c][t]) return 0.0;

    // HARD GUARD #1: same course already has a lecture at t → forbid
    if (slot_used[t * I->num_courses + c]) return 0.0;

    // HARD GUARD #2: clash with other courses already at t (curriculum/teacher)
    for (int cc = 0; cc < I->num_courses; ++cc) {
        if (slot_used[t * I->num_courses + cc] && I->conflict[c][cc])
            return 0.0;
    }

    // capacity scan among rooms FREE at time t
    int need = I->courses[c].students;
    int fit_free = 0;                     // any free room that FITS at t?
    int best_cap = 1000000000;            // smallest capacity among fitting free rooms
    int min_over = 1000000000;            // least overflow among free rooms (if none fit)

    for (int r = 0; r < I->num_rooms; ++r) {
        if (occ_room_time[r][t] > 0) continue;   // room already used at (r,t)
        int cap  = I->rooms[r].capacity;
        int over = need - cap;                    // <=0 means fits
        if (over <= 0) {
            fit_free = 1;
            if (cap < best_cap) best_cap = cap;  // prefer tight fits
        } else if (over < min_over) {
            min_over = over;
        }
    }

    // capacity-driven heuristic
    double eta;
    if (fit_free) {
        // prefer smaller fitting rooms
        eta = 6.0 + 3.0 / (1.0 + (double)best_cap);
    } else {
        // no free fitting room → strongly discourage
        double of = (double)(min_over > 0 ? min_over : 0);
        eta = 1e-9 / (1.0 + 0.03 * of);
    }
    if (eta < 1e-12) eta = 1e-12;

    return pow(tauLt, P.alpha) * pow(eta, P.beta);
}


// -------------------- main ACO kernel --------------------

// ---------- soft local search / post-optimizer (greedy, cheap) ----------

// overflow for (course c, room r)
static inline int overCR(const Instance* I, int c, int r){
    return I->courses[c].students - I->rooms[r].capacity; // <=0 is good
}

// recompute a lightweight room-time occupancy map from S
static void build_occ_rt(const Instance* I, const Timetable* S, int occ[][MAX_TIMESLOTS]){
    memset(occ, 0, sizeof(int)*MAX_ROOMS*MAX_TIMESLOTS);
    for(int L=0; L<I->total_lectures; ++L){
        int t = S->lec_to_timeslot[L]; if(t<0) continue;
        int r = S->lec_to_room[L];
        occ[r][t]++;
    }
}

// try improve capacity: (1) move to free room at same t, (2) swap rooms with another lecture at t
// try improve capacity: (1) move to free room at same t, (2) swap rooms with another lecture at t
static int improve_capacity_once(const Instance* I, Timetable* S){
    // TRACE("    soft: capacity start");
    const int R = I->num_rooms;
    const int T = I->num_timeslots;

    // heap occupancy: occ[r*T + t]
    int *occ = (int*)calloc((size_t)R * (size_t)T, sizeof(int));
    // if(!occ){ TRACE("    soft: capacity: alloc occ failed"); return 0; }

    // build occupancy from S
    for(int L=0; L<I->total_lectures; ++L){
        int t = S->lec_to_timeslot[L]; if(t<0 || t>=T) continue;
        int r = S->lec_to_room[L];     if(r<0 || r>=R) continue;
        occ[r*T + t] += 1;
    }

    int improved = 0;

    // (1) single reassign to a better free room
    for(int L=0; L<I->total_lectures; ++L){
        int t=S->lec_to_timeslot[L]; if(t<0 || t>=T) continue;
        int r=S->lec_to_room[L];     if(r<0 || r>=R) continue;
        int c=I->lecture_to_course[L];

        int best_r=r, best_over=I->courses[c].students - I->rooms[r].capacity;
        for(int rr=0; rr<R; ++rr){
            if(rr==r) continue;
            if(occ[rr*T + t] > 0) continue;                    // room busy at t
            int o = I->courses[c].students - I->rooms[rr].capacity; // <=0 fits
            // prefer any fit; otherwise least overflow (tie-break: smaller capacity)
            if( (o<=0 && best_over>0)
             || (o<=0 && best_over<=0 && I->rooms[rr].capacity < I->rooms[best_r].capacity)
             || (best_over>0 && o<best_over) ){
                best_r=rr; best_over=o;
            }
        }
        if(best_r!=r){
            occ[r*T + t]--; occ[best_r*T + t]++;
            S->lec_to_room[L]=best_r;
            improved=1;
        }
    }

    // (2) pairwise swap rooms at same timeslot to reduce total overflow
    for(int t=0; t<T; ++t){
        for(int r1=0; r1<R; ++r1){
            if(occ[r1*T + t]==0) continue;
            for(int r2=r1+1; r2<R; ++r2){
                if(occ[r2*T + t]==0) continue;

                // find the lectures currently in (r1,t) and (r2,t)
                int L1=-1, L2=-1;
                for(int L=0; L<I->total_lectures && (L1==-1 || L2==-1); ++L){
                    if(S->lec_to_timeslot[L]!=t) continue;
                    if(S->lec_to_room[L]==r1 && L1==-1) L1=L;
                    else if(S->lec_to_room[L]==r2 && L2==-1) L2=L;
                }
                if(L1==-1 || L2==-1) continue;

                int c1=I->lecture_to_course[L1], c2=I->lecture_to_course[L2];
                int cur = 0, swp = 0;
                int o11 = I->courses[c1].students - I->rooms[r1].capacity; if(o11>0) cur += o11;
                int o22 = I->courses[c2].students - I->rooms[r2].capacity; if(o22>0) cur += o22;
                int o12 = I->courses[c1].students - I->rooms[r2].capacity; if(o12>0) swp += o12;
                int o21 = I->courses[c2].students - I->rooms[r1].capacity; if(o21>0) swp += o21;

                if(swp < cur){
                    // perform swap
                    S->lec_to_room[L1]=r2;
                    S->lec_to_room[L2]=r1;
                    improved=1;
                }
            }
        }
    }

    free(occ);
    // TRACE("    soft: capacity end (improved=%d)", improved);
    return improved;
}


// Move up to 'need-used' lectures of a course to new days, cheapest first.
static int improve_mwd_once(const Instance* I, Timetable* S){
    // TRACE("    soft: mwd start");
    int changed = 0;

    for(int c=0; c<I->num_courses; ++c){
        int start = I->course_first_lecture[c];
        int cnt   = I->course_lecture_count[c];
        int need  = I->courses[c].min_working_days;
        if (cnt <= 1 || need <= 1) continue;

        // current day usage
        int daymask = 0;
        for(int k=0;k<cnt;++k){
            int L = start+k, t = S->lec_to_timeslot[L];
            if(t>=0) daymask |= 1 << (t / I->periods_per_day);
        }
        int used=0; for(int d=0; d<I->days; ++d) if(daymask & (1<<d)) ++used;
        if(used >= need) continue;

        int to_move = need - used;                 // how many new days we need
        int budget  = (to_move > 2 ? 2 : to_move); // cap per course per pass (tuning knob)

        // try to move 'budget' distinct lectures to distinct new days
        while(budget-- > 0){
            int best_L=-1, best_tt=-1, best_rr=-1, best_over=1e9, best_newday=-1;

            // choose which lecture to move (pick the cheapest relocation that creates a new day)
            for(int k=0;k<cnt;++k){
                int L = start+k;
                int t0 = S->lec_to_timeslot[L]; if(t0<0) continue;
                for(int d=0; d<I->days; ++d){
                    if(daymask & (1<<d)) continue; // we want a NEW day
                    for(int p=0; p<I->periods_per_day; ++p){
                        int tt = d*I->periods_per_day + p;
                        if(!I->feasible[c][tt]) continue;

                        // avoid conflicts at tt
                        int clash=0;
                        for(int L2=0; L2<I->total_lectures; ++L2){
                            if(S->lec_to_timeslot[L2]==tt && I->conflict[c][ I->lecture_to_course[L2] ]){ clash=1; break; }
                        }
                        if(clash) continue;

                        // pick any room free at tt (least overflow)
                        int best_local_rr=-1, best_local_over=1e9;
                        for(int rr=0; rr<I->num_rooms; ++rr){
                            int taken=0;
                            for(int L2=0; L2<I->total_lectures; ++L2)
                                if(S->lec_to_timeslot[L2]==tt && S->lec_to_room[L2]==rr){ taken=1; break; }
                            if(taken) continue;
                            int over = I->courses[c].students - I->rooms[rr].capacity;
                            if(over < best_local_over){ best_local_over=over; best_local_rr=rr; }
                        }
                        if(best_local_rr<0) continue;

                        // global best among all candidates
                        if(best_local_over < best_over){
                            best_over = best_local_over;
                            best_rr   = best_local_rr;
                            best_tt   = tt;
                            best_L    = L;
                            best_newday = d;
                        }
                    }
                }
            }

            if(best_L >= 0){
                S->lec_to_timeslot[best_L] = best_tt;
                S->lec_to_room[best_L]     = best_rr;
                daymask |= 1 << best_newday;
                changed = 1;
            }else{
                break; // no feasible new-day move found
            }
        }
    }
    // TRACE("    soft: mwd end (changed=%d)", changed);
    return changed;
}


// reduce Curriculum Compactness: move isolated singletons to neighbor periods
static int improve_cc_once(const Instance* I, Timetable* S){
    // TRACE("    soft: cc start (num_curricula=%d)", I->num_curricula);
    int changed=0;
    for(int cur=0; cur<I->num_curricula; ++cur){
        for(int d=0; d<I->days; ++d){
            for(int p=0; p<I->periods_per_day; ++p){
                int t=d*I->periods_per_day+p;
                // find lectures of this curriculum at day d
                // detect a singleton at (d,p): no lecture at p-1 or p+1 but one at p
                int has_prev=0, has_next=0, lone=0, Lsingleton=-1;
                for(int L=0; L<I->total_lectures; ++L){
                    int c=I->lecture_to_course[L];
                    int tt=S->lec_to_timeslot[L];
                    if(tt<0) continue;
                    if(tt==t) { lone=1; Lsingleton=L; }
                    if(p>0 && tt==t-1) has_prev=1;
                    if(p+1<I->periods_per_day && tt==t+1) has_next=1;
                }
                if(lone && !has_prev && !has_next){
                    // try move to p-1 or p+1 if feasible and non-conflicting
                    int c = I->lecture_to_course[Lsingleton];
                    int candidates[2]={ (p>0)? t-1 : -1, (p+1<I->periods_per_day)? t+1 : -1 };
                    for(int i=0;i<2;++i){
                        int tt=candidates[i]; if(tt<0) continue;
                        if(!I->feasible[c][tt]) continue;
                        // avoid conflicts at tt; find any free room (least overflow)
                        int clash=0;
                        for(int L2=0; L2<I->total_lectures; ++L2)
                            if(S->lec_to_timeslot[L2]==tt && I->conflict[c][ I->lecture_to_course[L2] ]){ clash=1; break; }
                        if(clash) continue;
                        int best_rr=-1, best_over=1e9;
                        for(int rr=0; rr<I->num_rooms; ++rr){
                            int taken=0;
                            for(int L2=0; L2<I->total_lectures; ++L2)
                                if(S->lec_to_timeslot[L2]==tt && S->lec_to_room[L2]==rr){ taken=1; break; }
                            if(taken) continue;
                            int o=overCR(I,c,rr);
                            if(o<best_over){ best_over=o; best_rr=rr; }
                        }
                        if(best_rr>=0){
                            S->lec_to_timeslot[Lsingleton]=tt;
                            S->lec_to_room[Lsingleton]=best_rr;
                            changed=1; break;
                        }
                    }
                }
            }
        }
    }
    // TRACE("    soft: cc end (changed=%d)", changed);
    return changed;
}

// prefer room stability for each course (assign same room when capacity allows)
static int improve_room_stability_once(const Instance* I, Timetable* S){
    // TRACE("    soft: rs start");
    int changed=0;
    for(int c=0;c<I->num_courses;++c){
        // pick a "home" room = smallest that fits (by capacity)
        int home=-1, best_cap=1e9;
        for(int r=0;r<I->num_rooms;++r){
            int o=overCR(I,c,r);
            if(o<=0 && I->rooms[r].capacity<best_cap){ best_cap=I->rooms[r].capacity; home=r; }
        }
        if(home<0) continue;
        // try to move lectures of c into 'home' if room is free at those times
        int start=I->course_first_lecture[c], cnt=I->course_lecture_count[c];
        for(int k=0;k<cnt;++k){
            int L=start+k, t=S->lec_to_timeslot[L]; if(t<0) continue;
            if(S->lec_to_room[L]==home) continue;
            int taken=0;
            for(int L2=0; L2<I->total_lectures; ++L2)
                if(L2!=L && S->lec_to_timeslot[L2]==t && S->lec_to_room[L2]==home){ taken=1; break; }
            if(!taken){
                S->lec_to_room[L]=home; changed=1;
            }
        }
    }
    // TRACE("    soft: rs end (changed=%d)", changed);
    return changed;
}

// full soft improvement (few cheap passes)
static void soft_post_opt(const Instance* I, Timetable* S){
    // TRACE("    soft: pass loop start");
    for(int it=0; it<SOFT_PASSES; ++it){
        // TRACE("    soft: pass=%d", it);
        int b = improve_mwd_once(I,S);      // spread days first
        int c = SOFT_ENABLE_CC ? improve_cc_once(I,S) : 0;
        int a = improve_capacity_once(I,S); // then fix room overflows
        int d = improve_room_stability_once(I,S);
        
        if(!(a||b||c||d)) break;
    }
    // TRACE("    soft: pass loop end");
}

// Fix any remaining hard violations: infeasible slots, clashes, room collisions.
static int hard_sanitize(const Instance* I, Timetable* S){
    int fixes = 0;

    // 1) move any lecture scheduled in an infeasible timeslot
    for(int L=0; L<I->total_lectures; ++L){
        int t = S->lec_to_timeslot[L];
        if(t<0) continue;
        int c = I->lecture_to_course[L];
        if(!I->feasible[c][t]){
            // move to first feasible slot with a free room (least overflow)
            int best_tt=-1, best_rr=-1, best_over=1e9;
            for(int tt=0; tt<I->num_timeslots; ++tt){
                if(!I->feasible[c][tt]) continue;
                // avoid teacher/curriculum clashes at tt
                int clash=0;
                for(int L2=0; L2<I->total_lectures; ++L2){
                    if(L2==L) continue;
                    if(S->lec_to_timeslot[L2]==tt){
                        int c2 = I->lecture_to_course[L2];
                        if(c==c2 || I->conflict[c][c2]){ clash=1; break; }
                    }
                }
                if(clash) continue;
                // pick any free room at tt (least overflow)
                for(int rr=0; rr<I->num_rooms; ++rr){
                    int taken=0;
                    for(int L2=0; L2<I->total_lectures; ++L2)
                        if(S->lec_to_timeslot[L2]==tt && S->lec_to_room[L2]==rr){ taken=1; break; }
                    if(taken) continue;
                    int over = I->courses[c].students - I->rooms[rr].capacity;
                    if(over < best_over){ best_over=over; best_tt=tt; best_rr=rr; }
                }
            }
            if(best_tt>=0){
                S->lec_to_timeslot[L]=best_tt;
                S->lec_to_room[L]=best_rr;
                ++fixes;
            }
        }
    }

    // 2) resolve room-time collisions (more than one lecture in same (room,t))
    //    try another free room at same t; else move to another feasible tt
    for(int t=0; t<I->num_timeslots; ++t){
        for(int r=0; r<I->num_rooms; ++r){
            // count occupants in (r,t)
            int occ_count=0, firstL=-1;
            for(int L=0; L<I->total_lectures; ++L){
                if(S->lec_to_timeslot[L]==t && S->lec_to_room[L]==r){
                    if(firstL<0) firstL=L;
                    ++occ_count;
                }
            }
            while(occ_count>1){ // move extra occupants
                // find one to move (not the first)
                int moveL=-1;
                for(int L=0; L<I->total_lectures; ++L)
                    if(S->lec_to_timeslot[L]==t && S->lec_to_room[L]==r && L!=firstL){ moveL=L; break; }
                if(moveL<0) break;

                int c = I->lecture_to_course[moveL];
                // try same t, different free room
                int moved=0;
                for(int rr=0; rr<I->num_rooms; ++rr){
                    if(rr==r) continue;
                    int taken=0;
                    for(int L2=0; L2<I->total_lectures; ++L2)
                        if(S->lec_to_timeslot[L2]==t && S->lec_to_room[L2]==rr){ taken=1; break; }
                    if(!taken){
                        S->lec_to_room[moveL]=rr; moved=1; ++fixes; --occ_count; break;
                    }
                }
                if(!moved){
                    // move to another feasible tt without clashes and with a free room
                    int best_tt=-1, best_rr=-1, best_over=1e9;
                    for(int tt=0; tt<I->num_timeslots; ++tt){
                        if(!I->feasible[c][tt]) continue;
                        // avoid teacher/curriculum clashes
                        int clash=0;
                        for(int L2=0; L2<I->total_lectures; ++L2){
                            if(L2==moveL) continue;
                            if(S->lec_to_timeslot[L2]==tt){
                                int c2 = I->lecture_to_course[L2];
                                if(c==c2 || I->conflict[c][c2]){ clash=1; break; }
                            }
                        }
                        if(clash) continue;
                        // any free room, least overflow
                        for(int rr=0; rr<I->num_rooms; ++rr){
                            int taken=0;
                            for(int L2=0; L2<I->total_lectures; ++L2)
                                if(S->lec_to_timeslot[L2]==tt && S->lec_to_room[L2]==rr){ taken=1; break; }
                            if(taken) continue;
                            int over = I->courses[c].students - I->rooms[rr].capacity;
                            if(over < best_over){ best_over=over; best_tt=tt; best_rr=rr; }
                        }
                    }
                    if(best_tt>=0){
                        S->lec_to_timeslot[moveL]=best_tt;
                        S->lec_to_room[moveL]=best_rr;
                        --occ_count; ++fixes;
                    }else{
                        break; // nothing we can do right now
                    }
                }
            }
        }
    }

    // 3) remove timeslot clashes (teacher/curriculum OR same course twice at same t)
    for(int t=0; t<I->num_timeslots; ++t){
        // mark a course already present at t
        // and move any extra lecture of same course to another slot
        for(int L=0; L<I->total_lectures; ++L){
            if(S->lec_to_timeslot[L]!=t) continue;
            int c = I->lecture_to_course[L];

            // if any OTHER lecture at t conflicts with course c → move this L
            int need_move = 0;
            for(int L2=0; L2<I->total_lectures; ++L2){
                if(L2==L) continue;
                if(S->lec_to_timeslot[L2]==t){
                    int c2 = I->lecture_to_course[L2];
                    if(c==c2 || I->conflict[c][c2]){ need_move=1; break; }
                }
            }
            if(!need_move) continue;

            // find another feasible tt without clashes (least overflow)
            int best_tt=-1, best_rr=-1, best_over=1e9;
            for(int tt=0; tt<I->num_timeslots; ++tt){
                if(tt==t) continue;
                if(!I->feasible[c][tt]) continue;
                int clash=0;
                for(int L3=0; L3<I->total_lectures; ++L3){
                    if(S->lec_to_timeslot[L3]==tt){
                        int c3 = I->lecture_to_course[L3];
                        if(c==c3 || I->conflict[c][c3]){ clash=1; break; }
                    }
                }
                if(clash) continue;
                for(int rr=0; rr<I->num_rooms; ++rr){
                    int taken=0;
                    for(int L3=0; L3<I->total_lectures; ++L3)
                        if(S->lec_to_timeslot[L3]==tt && S->lec_to_room[L3]==rr){ taken=1; break; }
                    if(taken) continue;
                    int over = I->courses[c].students - I->rooms[rr].capacity;
                    if(over < best_over){ best_over=over; best_tt=tt; best_rr=rr; }
                }
            }
            if(best_tt>=0){
                S->lec_to_timeslot[L]=best_tt;
                S->lec_to_room[L]=best_rr;
                ++fixes;
            }
        }
    }

    return fixes;
}


void aco_refine(Instance* I, ACOParams P, int seeds, const Timetable* seed_solutions,
                Timetable* out_best)
{
    // TRACE("enter ants=%d iters=%d (lectures=%d, timeslots=%d, rooms=%d)",
    //       P.ants, P.iters, I->total_lectures, I->num_timeslots, I->num_rooms);

    // --- allocate pheromone
    double** tau = new_matrix(I->total_lectures, I->num_timeslots, P.init_tau);
    if(!tau){
        // TRACE("alloc tau failed");
        *out_best = seed_solutions[0];
        return;
    }

    // --- seed pheromone from GA elites
    for(int s=0; s<seeds; ++s){
        const Timetable* S = &seed_solutions[s];
        for(int L=0; L<I->total_lectures; ++L){
            int t = S->lec_to_timeslot[L];
            if(t >= 0 && t < I->num_timeslots) tau[L][t] += 1.0;
        }
    }

    Timetable global_best = seed_solutions[0];

    // working buffers
    int*    cand  = (int*)malloc(sizeof(int) * I->num_timeslots);
    int*    order = (int*)malloc(sizeof(int) * I->total_lectures);
    double* diff  = (double*)malloc(sizeof(double) * I->total_lectures);
    if(!cand || !order || !diff){
        // TRACE("alloc work buffers failed (cand=%p order=%p diff=%p)", (void*)cand,(void*)order,(void*)diff);
        if(cand)  free(cand);
        if(order) free(order);
        if(diff)  free(diff);
        free_matrix(tau, I->total_lectures);
        *out_best = global_best;
        return;
    }

    // --- main loop
    for(int it=0; it<P.iters; ++it){
        // TRACE("iter=%d begin", it);

        // sort lectures by difficulty (hardest first)
        compute_difficulty(I, diff);
        make_order_desc(I, order, diff);

        for(int k_ant=0; k_ant<P.ants; ++k_ant){
            // if((k_ant % 16) == 0) TRACE("  iter=%d ant=%d", it, k_ant);

            Timetable T; memset(&T, 0, sizeof(T));
            for(int L=0; L<I->total_lectures; ++L){ T.lec_to_timeslot[L] = -1; T.lec_to_room[L] = 0; }

            // per-ant occupancy (heap to avoid big stacks)
            int (*occ_room_time)[MAX_TIMESLOTS] = calloc(I->num_rooms, sizeof *occ_room_time);
            if(!occ_room_time){
                // TRACE("    iter=%d ant=%d: alloc occ_room_time failed", it, k_ant);
                continue;
            }
            int *slot_used = calloc(I->num_timeslots * I->num_courses, sizeof(int)); // (t,course) map
            if(!slot_used){
                // TRACE("    iter=%d ant=%d: alloc slot_used failed", it, k_ant);
                free(occ_room_time);
                continue;
            }

            // construct solution in difficulty order
            for(int idx=0; idx<I->total_lectures; ++idx){
                int L = order[idx];

                // selective candidate set (top p% of feasible times)
                int m = build_selective_candidates(I, L, P.selective_p, cand, I->num_timeslots);
                if(m == 0){ T.lec_to_timeslot[L] = -1; continue; }
                if(m > 2048) m = 2048;  // guard scratch array

                // roulette weights with current occupancy/clash map
                double acc[2048], sum = 0.0;
                for(int i=0; i<m; ++i){
                    int t = cand[i];
                    double w = weight_for(P, I, L, t, tau[L][t], occ_room_time, slot_used);
                    sum += w; acc[i] = sum;
                }

                int pick = 0;
                if(sum > 0.0){
                    double r = unif01() * sum;
                    while(pick < m && r > acc[pick]) ++pick;
                    if(pick >= m) pick = m - 1;
                }else{
                    // fallback: first feasible
                    for(pick=0; pick<m; ++pick)
                        if(I->feasible[I->lecture_to_course[L]][cand[pick]]) break;
                    if(pick >= m) pick = m - 1;
                }

                int chosen_t = cand[pick];
                int c = I->lecture_to_course[L];

                int chosen_r = choose_room(I, c, chosen_t, occ_room_time);

                // commit
                T.lec_to_timeslot[L] = chosen_t;
                T.lec_to_room[L]     = chosen_r;
                // bounds guard (paranoid; indices should already be valid)
                if(chosen_r >= 0 && chosen_r < I->num_rooms &&
                   chosen_t >= 0 && chosen_t < I->num_timeslots){
                    occ_room_time[chosen_r][chosen_t]++;
                }
                slot_used[chosen_t * I->num_courses + c] = 1;

                // local pheromone update
                tau[L][chosen_t] = (1.0 - P.rho) * tau[L][chosen_t] + P.rho * 1.0;
            }

            // quick repairs (2 passes are usually enough to stay fast)
            repair_basic(I, &T);
            repair_basic(I, &T);
            
            // TRACE("  iter=%d ant=%d: soft_post_opt begin", it, k_ant);
            soft_post_opt(I, &T);
            // TRACE("  iter=%d ant=%d: soft_post_opt end", it, k_ant);
            
            for(int s=0; s<3; ++s){
                int f = hard_sanitize(I, &T);
                if(!f) break;
            }

            // score & keep global best
            // TRACE("  iter=%d ant=%d: evaluate", it, k_ant);
            evaluate(I, &T);
            if(dominates(&T, &global_best)) global_best = T;

            free(slot_used);
            free(occ_room_time);
        }

        // light evaporation (avoid stagnation)
        // TRACE("iter=%d evaporate", it);
        for(int L=0; L<I->total_lectures; ++L)
            for(int t=0; t<I->num_timeslots; ++t)
                tau[L][t] *= (1.0 - 0.01);

        // TRACE("iter=%d end", it);
    }

    TRACE("exit");
    *out_best = global_best;

    free(diff);
    free(order);
    free(cand);
    free_matrix(tau, I->total_lectures);
}
