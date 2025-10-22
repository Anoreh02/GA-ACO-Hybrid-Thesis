#pragma once
#include <stdint.h>

#define MAX_COURSES      512
#define MAX_ROOMS        128
#define MAX_INSTRUCTORS  256
#define MAX_DAYS         10
#define MAX_PERIODS      20
#define MAX_TIMESLOTS    (MAX_DAYS*MAX_PERIODS)   // timeslots = days * periods_per_day

#define MAX_CURRICULA        256
#define MAX_CURRIC_COURSES   512

// Upper bound for total lectures (safe for CB-CTT comp instances)
#define MAX_TOTAL_LECTURES   (MAX_COURSES * 16)

typedef struct {
    int id;                   // 0..num_courses-1
    int instructor_id;        // 0..num_instructors-1
    int room_type;            // CB-CTT: keep 0
    int hours_per_week;       // map from Lectures
    int is_lab;               // not used in CB-CTT
    int batch_id;             // not used in CB-CTT
    int min_working_days;     // from MinWorkingDays
    int students;             // Students
} Course;

typedef struct {
    int id;                   // 0..num_rooms-1
    int capacity;
    int room_type;            // CB-CTT: 0
} Room;

typedef struct {
    int id;                   // 0..num_instructors-1
    int pref[MAX_TIMESLOTS];  // CB-CTT: all zeros
    int target_load_hours;    // optional
} Instructor;

typedef struct {
    int days, periods_per_day;
    int num_courses, num_rooms, num_instructors, num_timeslots;

    Course     courses[MAX_COURSES];
    Room       rooms[MAX_ROOMS];
    Instructor instructors[MAX_INSTRUCTORS];

    // Feasibility for course c at timeslot t (accounts for UNAVAIL)
    uint8_t feasible[MAX_COURSES][MAX_TIMESLOTS];

    // Curricula
    int num_curricula;
    int curr_size[MAX_CURRICULA];
    int curr_courses[MAX_CURRICULA][MAX_CURRIC_COURSES];

    // Conflict matrix: 1 if courses conflict (same teacher OR same curriculum)
    uint8_t conflict[MAX_COURSES][MAX_COURSES];

    // Lecture indexing
    int total_lectures;
    int course_first_lecture[MAX_COURSES];   // offset in [0..total_lectures)
    int course_lecture_count[MAX_COURSES];   // lectures per course
    int lecture_to_course[MAX_TOTAL_LECTURES];

    // Optional summary
    int rooms_by_type[8][MAX_TIMESLOTS];

    // (optional) name maps (useful while parsing)
    int  num_course_ids, num_room_ids, num_teacher_ids;
    char course_name[MAX_COURSES][64];
    char room_name[MAX_ROOMS][64];
    char teacher_name[MAX_INSTRUCTORS][64];
} Instance;

// A timetable assigns each LECTURE to (timeslot, room)
typedef struct {
    int lec_to_timeslot[MAX_TOTAL_LECTURES]; // -1 if unassigned
    int lec_to_room   [MAX_TOTAL_LECTURES];  // room id

    // caches / scores
    int hard_violations;
    int soft_room_capacity;
    int soft_min_working_days;
    int soft_curriculum_compactness;
    int soft_room_stability;

    // extra metrics if you want to keep your tri-objective too
    double f1_conflicts;     // minimize (derived from hard/soft if desired)
    double f2_utilization;   // maximize
    double f3_workload_var;  // minimize
} Timetable;
