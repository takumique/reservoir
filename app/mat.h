#ifndef APP_MAT_H_
#define APP_MAT_H_

typedef float f32_t;
typedef double f64_t;

typedef struct {
  void *(*memory_alloc)(unsigned);
  void (*memory_free)(void *);
} mat_memory_t;

typedef struct {
  f32_t *data;
  unsigned n;
  unsigned m;
  unsigned t;
} mat_f32_t;

typedef struct {
  f64_t *data;
  unsigned n;
  unsigned m;
  unsigned t;
} mat_f64_t;

#define T(A) \
do { \
  (A).t = !(A).t; \
  unsigned _n = (A).n; \
  (A).n = (A).m; \
  (A).m = _n; \
} while(0)

#define _MAT(A, N, M) ((A).data + (A).m * N + M)
#define _MAT_T(A, N, M) ((A).data + (A).n * M + N)
#define MAT(A, N, M) ((A).t ? _MAT_T(A, N, M) : _MAT(A, N, M))

int mat_f32_new(mat_memory_t *mem, mat_f32_t *a, unsigned n, unsigned m);
void mat_f32_destroy(mat_memory_t *mem, mat_f32_t *a);
int mat_f32_copy(mat_f32_t *dst, mat_f32_t *src);
int mat_f32_zeros(mat_f32_t *a);
int mat_f32_sum(mat_f32_t *c, mat_f32_t *a, mat_f32_t *b);
int mat_f32_add(mat_f32_t *c, mat_f32_t *a, float l);
int mat_f32_product(mat_f32_t *c, mat_f32_t *a, mat_f32_t *b);
int mat_f32_mul(mat_f32_t *c, mat_f32_t *a, float l);
int mat_f32_identity(mat_f32_t *c, float l);
int mat_f32_inv(mat_f32_t *inv_a, mat_f32_t *a);
f32_t f32_random_normal(float mu, float sigma);
void mat_f32_random_normal(mat_f32_t *c, float mu, float sigma);
float mat_f32_max_abs_eigenval(mat_f32_t *a, mat_f32_t *x, mat_f32_t *y, unsigned lim);

int mat_f64_new(mat_memory_t *mem, mat_f64_t *a, unsigned n, unsigned m);
void mat_f64_destroy(mat_memory_t *mem, mat_f64_t *a);
int mat_f64_copy(mat_f64_t *dst, mat_f64_t *src);
int mat_f64_zeros(mat_f64_t *a);
int mat_f64_sum(mat_f64_t *c, mat_f64_t *a, mat_f64_t *b);
int mat_f64_add(mat_f64_t *c, mat_f64_t *a, double l);
int mat_f64_product(mat_f64_t *c, mat_f64_t *a, mat_f64_t *b);
int mat_f64_mul(mat_f64_t *c, mat_f64_t *a, double l);
int mat_f64_identity(mat_f64_t *c, double l);
int mat_f64_inv(mat_f64_t *inv_a, mat_f64_t *a);
f64_t f64_random_normal(double mu, double sigma);
void mat_f64_random_normal(mat_f64_t *c, double mu, double sigma);
double mat_f64_max_abs_eigenval(mat_f64_t *a, mat_f64_t *x, mat_f64_t *y, unsigned lim);

#if defined(PRECISION_F32)
#define MAT_NEW(...) mat_f32_new(__VA_ARGS__)
#define MAT_DESTROY(...) mat_f32_destroy(__VA_ARGS__)
#define MAT_COPY(...) mat_f32_copy(__VA_ARGS__)
#define MAT_ZEROS(...) mat_f32_zeros(__VA_ARGS__)
#define MAT_SUM(...) mat_f32_sum(__VA_ARGS__)
#define MAT_ADD(...) mat_f32_add(__VA_ARGS__)
#define MAT_PRODUCT(...) mat_f32_product(__VA_ARGS__)
#define MAT_MUL(...) mat_f32_mul(__VA_ARGS__)
#define MAT_IDENTITY(...) mat_f32_identity(__VA_ARGS__)
#define MAT_INV(...) mat_f32_inv(__VA_ARGS__)
#define RANDOM_NORMAL(...) f32_random_normal(__VA_ARGS__)
#define MAT_RANDOM_NORMAL(...) mat_f32_random_normal(__VA_ARGS__)
#define MAT_MAX_ABS_EIGENVAL(...) mat_f32_max_abs_eigenval(__VA_ARGS__)
#elif defined(PRECISION_F32)
#define MAT_NEW(...) mat_f64_new(__VA_ARGS__)
#define MAT_DESTROY(...) mat_f64_destroy(__VA_ARGS__)
#define MAT_ZEROS(...) mat_f64_zeros(__VA_ARGS__)
#define MAT_SUM(...) mat_f64_sum(__VA_ARGS__)
#define MAT_ADD(...) mat_f64_add(__VA_ARGS__)
#define MAT_PRODUCT(...) mat_f64_product(__VA_ARGS__)
#define MAT_MUL(...) mat_f64_mul(__VA_ARGS__)
#define MAT_IDENTITY(...) mat_f64_identity(__VA_ARGS__)
#define MAT_INV(...) mat_f64_inv(__VA_ARGS__)
#define RANDOM_NORMAL(...) f64_random_normal(__VA_ARGS__)
#define MAT_RANDOM_NORMAL(...) mat_f64_random_normal(__VA_ARGS__)
#define MAT_MAX_ABS_EIGENVAL(...) mat_f64_max_abs_eigenval(__VA_ARGS__)
#endif

#endif /* APP_MAT_H_ */
