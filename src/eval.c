#include <string.h>
#include <math.h>
#include "eval.h"
#include "util.h"

// helpers
static inline int day_of(const Instance* I, int t){ return t / I->periods_per_day; }
static inline int per_of(const Instance* I, int t){ return t % I->periods_per_day; }

static double workload_variance(Instance* I, const Timetable* S){
    int maxI = I->num_instructors;
    int hours[MAX_INSTRUCTORS]; memset(hours,0,sizeof(hours));
    for(int L=0; L<I->total_lectures; ++L){
        int t = S->lec_to_timeslot[L];
        if(t>=0){
            int c = I->lecture_to_course[L];
            hours[I->courses[c].instructor_id] += 1; // per-lecture = 1 period
        }
    }
    double mean=0; for(int i=0;i<maxI;++i) mean += hours[i]; mean/= (double)maxI;
    double var=0;  for(int i=0;i<maxI;++i){ double d=hours[i]-mean; var+=d*d; }
    return (maxI>0)? var/(double)maxI : 0.0;
}

static void compute_soft_penalties(Instance* I, Timetable* S){
    // Reset softs
    S->soft_room_capacity = 0;
    S->soft_min_working_days = 0;
    S->soft_curriculum_compactness = 0;
    S->soft_room_stability = 0;

    // RoomCapacity: sum over-occupations
    for(int L=0; L<I->total_lectures; ++L){
        int t = S->lec_to_timeslot[L];
        if(t<0) continue;
        int r = S->lec_to_room[L];
        int c = I->lecture_to_course[L];
        int over = I->courses[c].students - I->rooms[r].capacity;
        if(over > 0) S->soft_room_capacity += over;
    }

    // MinimumWorkingDays: 5 * shortfall
    for(int c=0;c<I->num_courses;++c){
        int start = I->course_first_lecture[c];
        int cnt   = I->course_lecture_count[c];
        int used_days_mask = 0;
        for(int k=0;k<cnt;++k){
            int L = start+k;
            int t = S->lec_to_timeslot[L];
            if(t>=0){
                int d = day_of(I,t);
                used_days_mask |= (1<<d);
            }
        }
        int used_days=0;
        for(int d=0; d<I->days; ++d) if(used_days_mask & (1<<d)) ++used_days;
        int shortfall = I->courses[c].min_working_days - used_days;
        if(shortfall>0) S->soft_min_working_days += 5*shortfall;
    }

    // RoomStability: (distinct rooms used per course) - 1
    for(int c=0;c<I->num_courses;++c){
        int used_room[ MAX_ROOMS ]; memset(used_room,0,sizeof(used_room));
        int start = I->course_first_lecture[c];
        int cnt   = I->course_lecture_count[c];
        for(int k=0;k<cnt;++k){
            int L = start+k;
            int t = S->lec_to_timeslot[L];
            if(t>=0){
                int r = S->lec_to_room[L];
                if(r>=0 && r<I->num_rooms) used_room[r]=1;
            }
        }
        int rooms_dist=0; for(int r=0;r<I->num_rooms;++r) rooms_dist += used_room[r];
        if(rooms_dist>1) S->soft_room_stability += (rooms_dist - 1);
    }

    // CurriculumCompactness: for each curriculum, day; penalize isolated lectures (no neighbor in adjacent periods)
    for(int cur=0; cur<I->num_curricula; ++cur){
        for(int d=0; d<I->days; ++d){
            // mark occupied periods for this curriculum and day
            int occ[MAX_PERIODS]; memset(occ,0,sizeof(occ));
            for(int idx=0; idx<I->curr_size[cur]; ++idx){
                int c = I->curr_courses[cur][idx];
                int start = I->course_first_lecture[c];
                int cnt   = I->course_lecture_count[c];
                for(int k=0;k<cnt;++k){
                    int L = start+k;
                    int t = S->lec_to_timeslot[L];
                    if(t>=0 && day_of(I,t)==d){
                        int p = per_of(I,t);
                        occ[p] = 1;
                    }
                }
            }
            // count isolated periods in this day
            for(int p=0; p<I->periods_per_day; ++p){
                if(!occ[p]) continue;
                int left  = (p>0)? occ[p-1] : 0;
                int right = (p+1<I->periods_per_day)? occ[p+1] : 0;
                if(!left && !right) S->soft_curriculum_compactness += 2;
            }
        }
    }
}

void evaluate(Instance* I, Timetable* S){
    // reset hard violations
    S->hard_violations = 0;

    // 1) Lectures: all lectures must be scheduled
    for(int L=0; L<I->total_lectures; ++L){
        if(S->lec_to_timeslot[L] < 0) S->hard_violations++;
    }

    // 2) RoomOccupancy: same room, same timeslot
    // build occupancy table [room][timeslot]
    static int occ[MAX_ROOMS][MAX_TIMESLOTS];
    memset(occ,0,sizeof(occ));
    for(int L=0; L<I->total_lectures; ++L){
        int t = S->lec_to_timeslot[L];
        if(t<0) continue;
        int r = S->lec_to_room[L];
        if(++occ[r][t] > 1) S->hard_violations++;
    }

    // 3) Conflicts: lectures of conflicting courses may not share same timeslot
    // (we don’t care about rooms here—any room counts)
    static int slot_courses[MAX_TIMESLOTS][MAX_COURSES];
    static int slot_count[MAX_TIMESLOTS];
    memset(slot_count,0,sizeof(slot_count));
    for(int t=0; t<I->num_timeslots; ++t) memset(slot_courses[t], -1, sizeof(int)*MAX_COURSES);

    for(int L=0; L<I->total_lectures; ++L){
        int t = S->lec_to_timeslot[L];
        if(t<0) continue;
        int c = I->lecture_to_course[L];
        slot_courses[t][ slot_count[t]++ ] = c;
    }
    for(int t=0; t<I->num_timeslots; ++t){
        for(int i=0;i<slot_count[t];++i){
            int a = slot_courses[t][i]; if(a<0) continue;
            for(int j=i+1;j<slot_count[t];++j){
                int b = slot_courses[t][j]; if(b<0) continue;
                if(I->conflict[a][b]) S->hard_violations++;
            }
        }
    }

    // 4) Availabilities: if infeasible[c][t]==0 but scheduled => violation
    for(int L=0; L<I->total_lectures; ++L){
        int t = S->lec_to_timeslot[L]; if(t<0) continue;
        int c = I->lecture_to_course[L];
        if(!I->feasible[c][t]) S->hard_violations++;
    }

    // Soft penalties
    compute_soft_penalties(I,S);

    // Optional diagnostic objectives
    // f1: you can combine hard + soft to a proxy; here keep as soft sum
    S->f1_conflicts = (double)ctt_soft_score(S);
    // f2 utilization: number of lectures placed
    int placed=0; for(int L=0; L<I->total_lectures; ++L) if(S->lec_to_timeslot[L]>=0) placed++;
    S->f2_utilization = (double)placed;
    // f3 workload variance
    S->f3_workload_var = workload_variance(I,S);
}

int dominates(const Timetable* A, const Timetable* B){
    // Multi-objective: minimize f1, f3; maximize f2
    int better_or_equal = (A->f1_conflicts <= B->f1_conflicts) &&
                          (A->f2_utilization >= B->f2_utilization) &&
                          (A->f3_workload_var <= B->f3_workload_var);
    int strictly_better = (A->f1_conflicts <  B->f1_conflicts) ||
                          (A->f2_utilization >  B->f2_utilization) ||
                          (A->f3_workload_var <  B->f3_workload_var);
    return better_or_equal && strictly_better;
}

int ctt_soft_score(const Timetable* S){
    return S->soft_room_capacity
         + S->soft_min_working_days
         + S->soft_curriculum_compactness
         + S->soft_room_stability;
}
