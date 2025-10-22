#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>
#include "ctt.h"

static void trim(char* s){
    char* p=s; while(isspace((unsigned char)*p)) p++;
    if(p!=s) memmove(s,p,strlen(p)+1);
    size_t n=strlen(s);
    while(n>0 && isspace((unsigned char)s[n-1])) s[--n]='\0';
}

static int course_index(Instance* I, const char* name){
    for(int i=0;i<I->num_course_ids;i++) if(strcmp(I->course_name[i], name)==0) return i;
    if(I->num_course_ids>=MAX_COURSES) return -1;
    strncpy(I->course_name[I->num_course_ids], name, 63);
    I->course_name[I->num_course_ids][63]='\0';
    return I->num_course_ids++;
}

static int room_index(Instance* I, const char* name){
    for(int i=0;i<I->num_room_ids;i++) if(strcmp(I->room_name[i], name)==0) return i;
    if(I->num_room_ids>=MAX_ROOMS) return -1;
    strncpy(I->room_name[I->num_room_ids], name, 63);
    I->room_name[I->num_room_ids][63]='\0';
    return I->num_room_ids++;
}

static int teacher_index(Instance* I, const char* name){
    for(int i=0;i<I->num_teacher_ids;i++) if(strcmp(I->teacher_name[i], name)==0) return i;
    if(I->num_teacher_ids>=MAX_INSTRUCTORS) return -1;
    strncpy(I->teacher_name[I->num_teacher_ids], name, 63);
    I->teacher_name[I->num_teacher_ids][63]='\0';
    return I->num_teacher_ids++;
}

int load_ctt(const char* path, Instance* I){
    FILE* f = fopen(path, "r");
    if(!f) return -1;

    memset(I, 0, sizeof(*I));

    char line[1024];
    int expected_courses=0, expected_rooms=0, expected_curricula=0, expected_constraints=0;

    enum Section { SEC_NONE, SEC_COURSES, SEC_ROOMS, SEC_CURRICULA, SEC_UNAVAIL } sec = SEC_NONE;

    while(fgets(line,sizeof(line),f)){
        trim(line);
        if(line[0]=='\0' || line[0]=='%') continue;

        if(strncmp(line,"Name:",5)==0) continue;
        if(sscanf(line,"Days: %d",&I->days)==1) continue;
        if(sscanf(line,"Periods_per_day: %d",&I->periods_per_day)==1) continue;
        if(sscanf(line,"Courses: %d",&expected_courses)==1) continue;
        if(sscanf(line,"Rooms: %d",&expected_rooms)==1) continue;
        if(sscanf(line,"Curricula: %d",&expected_curricula)==1) continue;
        if(sscanf(line,"Constraints: %d",&expected_constraints)==1) continue;

        if(strcmp(line,"COURSES:")==0){ sec=SEC_COURSES; continue; }
        if(strcmp(line,"ROOMS:")==0){ sec=SEC_ROOMS; continue; }
        if(strcmp(line,"CURRICULA:")==0){ sec=SEC_CURRICULA; continue; }
        if(strcmp(line,"UNAVAILABILITY_CONSTRAINTS:")==0){ sec=SEC_UNAVAIL; continue; }

        if(sec==SEC_COURSES){
            // CourseID TeacherID Lectures MinWorkingDays Students
            char cID[64], tID[64]; int lect, mindays, students;
            if(sscanf(line,"%63s %63s %d %d %d", cID, tID, &lect, &mindays, &students)==5){
                int c = course_index(I,cID);
                int inst = teacher_index(I,tID);
                if(c<0 || inst<0){ fclose(f); return -2; }

                I->courses[c].id = c;
                I->courses[c].instructor_id = inst;
                I->courses[c].hours_per_week = lect;
                I->courses[c].min_working_days = mindays;
                I->courses[c].students = students;
                I->courses[c].room_type = 0;
                I->num_courses = I->num_course_ids;
            }
            continue;
        }

        if(sec==SEC_ROOMS){
            // RoomID Capacity
            char rID[64]; int cap;
            if(sscanf(line,"%63s %d", rID, &cap)==2){
                int r = room_index(I,rID);
                if(r<0){ fclose(f); return -3; }
                I->rooms[r].id = r;
                I->rooms[r].capacity = cap;
                I->rooms[r].room_type = 0;
                I->num_rooms = I->num_room_ids;
            }
            continue;
        }

        if(sec==SEC_CURRICULA){
            // CurriculumID N course1 course2 ... courseN
            char tmp[1024]; strncpy(tmp,line,sizeof(tmp)-1); tmp[sizeof(tmp)-1]='\0';
            char* tok = strtok(tmp," \t");
            if(!tok) continue; // curID
            tok = strtok(NULL," \t"); if(!tok) continue; // N
            int n = atoi(tok);
            int idx = I->num_curricula++;
            I->curr_size[idx]=0;
            for(int i=0;i<n;i++){
                tok = strtok(NULL," \t");
                if(!tok) break;
                int c = course_index(I, tok);
                I->curr_courses[idx][ I->curr_size[idx]++ ] = c;
            }
            continue;
        }

        if(sec==SEC_UNAVAIL){
            // CourseID Day Period
            char cID[64]; int d,pd;
            if(sscanf(line,"%63s %d %d", cID, &d, &pd)==3){
                int c = course_index(I, cID);
                // num_timeslots finalization later
                // We'll mark feasibility after we know num_timeslots
                // For now, stash in a temp array? Simpler: set once timeslots known.
                // Here we assume Days/Periods_per_day already parsed.
                if(I->num_timeslots==0){
                    I->num_timeslots = I->days * I->periods_per_day;
                    for(int cc=0; cc<MAX_COURSES; ++cc)
                        for(int t=0;t<MAX_TIMESLOTS;++t)
                            I->feasible[cc][t] = (t<I->num_timeslots);
                }
                int t = d*I->periods_per_day + pd;
                if(c>=0 && t>=0 && t<I->num_timeslots) I->feasible[c][t]=0;
            }
            continue;
        }
    }
    fclose(f);

    // finalize counts
    if(I->num_timeslots==0){
        I->num_timeslots = I->days * I->periods_per_day;
        for(int c=0;c<I->num_courses;++c)
            for(int t=0;t<I->num_timeslots;++t)
                I->feasible[c][t]=1;
    }

    // instructors
    I->num_instructors = I->num_teacher_ids;
    for(int i=0;i<I->num_instructors;++i){
        I->instructors[i].id=i;
        I->instructors[i].target_load_hours=0;
        for(int t=0;t<I->num_timeslots;++t) I->instructors[i].pref[t]=0;
    }

    // rooms_by_type (all 0)
    for(int t=0;t<I->num_timeslots;++t) I->rooms_by_type[0][t] = I->num_rooms;

    // Build conflict matrix (teacher)
    memset(I->conflict,0,sizeof(I->conflict));
    for(int a=0; a<I->num_courses; ++a)
        for(int b=a+1; b<I->num_courses; ++b)
            if(I->courses[a].instructor_id == I->courses[b].instructor_id)
                I->conflict[a][b] = I->conflict[b][a] = 1;

    // Add curriculum-based conflicts
    for(int k=0; k<I->num_curricula; ++k){
        for(int i=0; i<I->curr_size[k]; ++i){
            int a = I->curr_courses[k][i];
            for(int j=i+1; j<I->curr_size[k]; ++j){
                int b = I->curr_courses[k][j];
                I->conflict[a][b] = I->conflict[b][a] = 1;
            }
        }
    }

    // Build lecture indexing
    I->total_lectures = 0;
    for(int c=0;c<I->num_courses;++c){
        int L = I->courses[c].hours_per_week;
        I->course_first_lecture[c]  = I->total_lectures;
        I->course_lecture_count[c]  = L;
        for(int k=0;k<L;++k){
            int gid = I->total_lectures + k;
            I->lecture_to_course[gid] = c;
        }
        I->total_lectures += L;
    }

    return 0;
}
