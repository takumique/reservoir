#ifndef APP_LINALG_H_
#define APP_LINALG_H_

void matrix_product_f32(float *c, float *a, float *b, unsigned n, unsigned m, unsigned p);
void matrix_product_left_t_f32(float *c, float *a, float *b, unsigned n, unsigned m, unsigned p);
void matrix_product_right_t_f32(float *c, float *a, float *b, unsigned n, unsigned m, unsigned p);
void vector_multiplication_f32(float *c, float *a, float l, unsigned n);
void matrix_multiplication_f32(float *c, float *a, float l, unsigned n, unsigned m);
void identity_multiplied_f32(float *c, float l, unsigned n);
void identity_f32(float *c, unsigned n);
void inv_f32(float *inv_a, float *a, unsigned n);
float random_normal_f32(float mu, float sigma);
void vector_random_normal_f32(float *a, unsigned n, float mu, float sigma);
void matrix_random_normal_f32(float *a, unsigned n, unsigned m, float mu, float sigma);
float max_abs_eigenval_f32(float *a, float *x, float *y, unsigned n, unsigned lim);

void matrix_product_f64(double *c, double *a, double *b, unsigned n, unsigned m, unsigned p);
void matrix_product_left_t_f64(double *c, double *a, double *b, unsigned n, unsigned m, unsigned p);
void matrix_product_right_t_f64(double *c, double *a, double *b, unsigned n, unsigned m, unsigned p);
void vector_multiplication_f64(double *c, double *a, double l, unsigned n);
void matrix_multiplication_f64(double *c, double *a, double l, unsigned n, unsigned m);
void identity_multiplied_f64(double *c, double l, unsigned n);
void identity_f64(double *c, unsigned n);
void inv_f64(double *inv_a, double *a, unsigned n);
double random_normal_f64(double mu, double sigma);
void vector_random_normal_f64(double *a, unsigned n, double mu, double sigma);
void matrix_random_normal_f64(double *a, unsigned n, unsigned m, double mu, double sigma);
double max_abs_eigenval_f64(double *a, double *x, double *y, unsigned n, unsigned lim);

#ifdef SINGLE_PRECISION
#define MATRIX_PRODUCT(...) matrix_product_f32(__VA_ARGS__)
#define MATRIX_PRODUCT_LEFT_T(...) matrix_product_left_t_f32(__VA_ARGS__)
#define MATRIX_PRODUCT_RIGHT_T(...) matrix_product_right_t_f32(__VA_ARGS__)
#define VECTOR_MULTIPLICATION(...) vector_multiplication_f32(__VA_ARGS__)
#define MATRIX_MULTIPLICATION(...) matrix_multiplication_f32(__VA_ARGS__)
#define IDENTITY_MULTIPLIED(...) identity_multiplied_f32(__VA_ARGS__)
#define IDENTITY(...) identity_f32(__VA_ARGS__)
#define INV(...) inv_f32(__VA_ARGS__)
#define RANDOM_NORMAL(...) random_normal_f32(__VA_ARGS__)
#define VECTOR_RANDOM_NORMAL(...) vector_random_normal_f32(__VA_ARGS__)
#define MATRIX_RANDOM_NORMAL(...) matrix_random_normal_f32(__VA_ARGS__)
#define MAX_ABS_EIGENVAL(...) max_abs_eigenval_f32(__VA_ARGS__)
#else
#define MATRIX_PRODUCT(...) matrix_product_f64(__VA_ARGS__)
#define MATRIX_PRODUCT_LEFT_T(...) matrix_product_left_t_f64(__VA_ARGS__)
#define MATRIX_PRODUCT_RIGHT_T(...) matrix_product_right_t_f64(__VA_ARGS__)
#define VECTOR_MULTIPLICATION(...) vector_multiplication_f64(__VA_ARGS__)
#define MATRIX_MULTIPLICATION(...) matrix_multiplication_f64(__VA_ARGS__)
#define IDENTITY_MULTIPLIED(...) identity_multiplied_f64(__VA_ARGS__)
#define IDENTITY(...) identity_f64(__VA_ARGS__)
#define INV(...) inv_f64(__VA_ARGS__)
#define RANDOM_NORMAL(...) random_normal_f64(__VA_ARGS__)
#define VECTOR_RANDOM_NORMAL(...) vector_random_normal_f64(__VA_ARGS__)
#define MATRIX_RANDOM_NORMAL(...) matrix_random_normal_f64(__VA_ARGS__)
#define MAX_ABS_EIGENVAL(...) max_abs_eigenval_f64(__VA_ARGS__)
#endif

#endif /* APP_LINALG_H_ */
