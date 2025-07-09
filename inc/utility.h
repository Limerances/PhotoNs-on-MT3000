#ifndef UTILITY_H
#define UTILITY_H

#ifndef LOGLEVEL
#define LOGLEVEL 0  // 0: only error 1: error and warning 2: error, warning and info 3: all
#endif


#define LOG_ERR 0
#define LOG_WRN 1
#define LOG_INF 2
#define LOG_DBG 3
 
#define LOG(lvl, ...) \
    do { \
        if(lvl <= LOGLEVEL) { \
            switch(lvl) { \
                case LOG_ERR: \
                              fprintf(stderr, "\"%s\" %d [err]: ", __FILE__, __LINE__); \
                break; \
                case LOG_WRN: \
                              fprintf(stderr, "\"%s\" %d [wrn]: ", __FILE__, __LINE__); \
                break; \
                case LOG_INF: \
                              fprintf(stderr, "\"%s\" %d [info]: ", __FILE__, __LINE__); \
                break; \
                case LOG_DBG: \
                              fprintf(stderr, "\"%s\" %d [dbg]: ", __FILE__, __LINE__); \
                break; \
                default: \
                         fprintf(stderr, "\"%s\" %d [???]: ", __FILE__, __LINE__); \
                break; \
            } \
            fprintf(stderr, __VA_ARGS__); \
            fprintf(stderr, "\n"); \
        } \
    } while(0)

#ifdef __cplusplus
extern "C" {
#endif
double dtime();
void make_dir(const char* path, const char* name, int r_snap);
#ifdef __cplusplus
}
#endif

#endif
