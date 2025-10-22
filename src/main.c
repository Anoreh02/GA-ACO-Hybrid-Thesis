#include <stdio.h>
#include <stdlib.h>
#include "model.h"
#include "eval.h"
#include "ga.h"
#include "aco.h"
#include "ctt.h"
#include "util.h"

enum { K = 8 };

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
    ACOParams ap = {.ants=40, .iters=120, .alpha=1.0, .beta=4.5, .rho=0.2, .init_tau=0.1, .selective_p=0.30};

    double t0 = now_sec();

    printf("[MAIN] Starting GA...\n"); fflush(stdout);
    Timetable best_ga;
    Timetable seeds[K];
    double t_ga0 = now_sec();
    ga_run(&I, gp, &best_ga, K, seeds);
    double t_ga1 = now_sec();

    printf("[MAIN] GA done. (%.3fs)\n", t_ga1 - t_ga0);

    printf("[MAIN] Starting ACO...\n"); fflush(stdout);
    Timetable best_final;
    double t_aco0 = now_sec();
    aco_refine(&I, ap, K, seeds, &best_final);
    double t_aco1 = now_sec();

    printf("[MAIN] ACO done. (%.3fs)\n", t_aco1 - t_aco0);


    printf("[TIME] total=%.3fs | GA=%.3fs | ACO=%.3fs\n",
       (t_aco1 - t0), (t_ga1 - t_ga0), (t_aco1 - t_aco0));
       
    evaluate(&I,&best_final);
    printf("Hard=%d | Soft=%d (Cap=%d, MWD=%d, CC=%d, RS=%d)\n",
           best_final.hard_violations,
           ctt_soft_score(&best_final),
           best_final.soft_room_capacity,
           best_final.soft_min_working_days,
           best_final.soft_curriculum_compactness,
           best_final.soft_room_stability);

    printf("f1=%.2f f2=%.2f f3=%.4f\n",
           best_final.f1_conflicts, best_final.f2_utilization, best_final.f3_workload_var);
    return 0;
}
