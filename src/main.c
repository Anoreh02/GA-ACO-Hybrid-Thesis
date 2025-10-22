#include <stdio.h>
#include <stdlib.h>
#include "model.h"
#include "eval.h"
#include "ga.h"
#include "aco.h"
#include "ctt.h"
#include "util.h"

enum { K = 8 };

void export_timetable_csv(const Instance* I, const Timetable* T, const char* filename){
    FILE* f = fopen(filename, "w");
    if(!f){ perror("fopen"); return; }
    fprintf(f, "Course,Instructor,Room,Day,Period\n");
    for(int L = 0; L < I->total_lectures; ++L){
        int c = I->lecture_to_course[L];
        int t = T->lec_to_timeslot[L];
        int r = T->lec_to_room[L];
        int day = t / I->periods_per_day;
        int period = t % I->periods_per_day;
        fprintf(f, "%d,%d,%d,%d,%d\n", c, I->courses[c].instructor_id, r, day+1, period+1);
    }
    fclose(f);
    printf("[MAIN] Timetable exported to %s\n", filename);
}

int main(int argc, char** argv){
    srand(42);

    const char* path = (argc>=2)? argv[1] : "data/comp01.ctt";

    Instance I;
    if(load_ctt(path, &I)!=0){
        fprintf(stderr,"[ERR] Failed to load CTT file: %s\n", path);
        return 1;
    }
    printf("Loaded %s | courses=%d rooms=%d instructors=%d days=%d periods=%d lectures=%d\n",
           path, I.num_courses, I.num_rooms, I.num_instructors, I.days, I.periods_per_day, I.total_lectures);
    fflush(stdout);

    // Small, fast settings first so you see output quickly. Tune up later.
    GAParams gp = {.pop_size=40, .iters=30, .pm=0.05};
    ACOParams ap = {.ants=32, .iters=100, .alpha=1.0, .beta=4.5, .rho=0.2, .init_tau=0.1, .selective_p=0.30};

    double t0 = now_sec();

    printf("[MAIN] Starting GA...\n"); fflush(stdout);
    Timetable best_ga;
    Timetable seeds[K];
    double t_ga0 = now_sec();
    ga_run(&I, gp, &best_ga, K, seeds);
    double t_ga1 = now_sec();

    printf("[MAIN] GA done. (%.3fs)\n", t_ga1 - t_ga0);

    // --- NEW: evaluate, print, and export GA best before ACO ---
    evaluate(&I, &best_ga);
    printf("[MAIN] GA best: f1=%.2f f2=%.2f f3=%.4f | Hard=%d | Soft=%d\n",
           best_ga.f1_conflicts, best_ga.f2_utilization, best_ga.f3_workload_var,
           best_ga.hard_violations, ctt_soft_score(&best_ga));
    export_timetable_csv(&I, &best_ga, "timetable_ga_seed.csv");
    // ------------------------------------------------------------

    printf("[MAIN] Starting ACO...\n"); fflush(stdout);
    Timetable best_final;
    double t_aco0 = now_sec();
    aco_refine(&I, ap, K, seeds, &best_final);
    double t_aco1 = now_sec();

    printf("[MAIN] ACO done. (%.3fs)\n", t_aco1 - t_aco0);

    // Evaluate final (needed for delta + final report)
    evaluate(&I, &best_final);

    // --- NEW: show ACO - GA deltas so you can see impact clearly ---
    printf("[DELTA] ACO - GA: f1=%.2f f2=%.2f f3=%.4f (negative f1/f3 is better, positive f2 is better)\n",
           best_final.f1_conflicts - best_ga.f1_conflicts,
           best_final.f2_utilization - best_ga.f2_utilization,
           best_final.f3_workload_var - best_ga.f3_workload_var);
    // ----------------------------------------------------------------

    printf("[TIME] total=%.3fs | GA=%.3fs | ACO=%.3fs\n",
       (t_aco1 - t0), (t_ga1 - t_ga0), (t_aco1 - t_aco0));

    printf("Hard=%d | Soft=%d (Cap=%d, MWD=%d, CC=%d, RS=%d)\n",
           best_final.hard_violations,
           ctt_soft_score(&best_final),
           best_final.soft_room_capacity,
           best_final.soft_min_working_days,
           best_final.soft_curriculum_compactness,
           best_final.soft_room_stability);

    printf("f1=%.2f f2=%.2f f3=%.4f\n",
           best_final.f1_conflicts, best_final.f2_utilization, best_final.f3_workload_var);

    export_timetable_csv(&I, &best_final, "timetable.csv");

    return 0;
}
