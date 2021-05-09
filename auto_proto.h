/* ab2centroid.c */
/* ab2scatter.c */
/* abadd.c */
/* abconvert.c */
void taperf(double *, int, int, double, int, int, double *, double *, double *, double);
void phelp(char *);
/* calc_radial_corr.c */
/* chebyshev.c */
double chebev(double, double, double *, int, double);
double norm_chebev(double, double, double *, int, double);
double chebnorm(int);
double cheb_norm_damping(int, int, int);
double chebev_der_int(double, double, int);
double norm_chebev_der_int(double, double, int);
/* cmodelcorr.c */
/* cmodellinreg.c */
/* cmodelmeancorr.c */
/* cmodelpower.c */
/* cmradialcorr.c */
/* convert_she_model.c */
/* cradialcorr.c */
/* determine_coeff.c */
void determine_coeff(double *, int, double *, double *, int, double, double, double, double, double *, double *, double *, double (*)(double, double, double *, int, double), double (*)(int, int, int));
void svd_solver(double *, double *, double *, int, int);
void svd_driver(double *, int, int, double *, double *);
/* dfour1.c */
void dfour1(double *, unsigned long, int);
/* extract_layer.c */
/* extract_model_depths.c */
/* fitxy.c */
void linreg_fit(double *, double *, int, double *, int, double *, double *, double *, double *, double *, double *);
void linreg_fitexy(double *, double *, int, double *, double *, double *, double *, double *, double *, double *, double *);
double chixy(double);
double brent(double, double, double, double (*)(void), double, double *);
double gammq(double, double);
void gcf(double *, double, double, double *);
void gser(double *, double, double, double *);
void avevar(double *, unsigned long, double *, double *);
void fit(double *, double *, int, double *, int, double *, double *, double *, double *, double *, double *);
void mnbrak(double *, double *, double *, double *, double *, double *, double (*)(void));
double zbrent(double (*)(void), double, double, double);
/* four1.c */
void four1(double *, unsigned long, int);
/* gauleg.c */
void gauleg_orig(double, double, double *, double *, int);
/* gsh_handling.c */
void write_gsh_coeff_set(int, int, double *, double *, double *, double *, double *, double *, double, FILE *);
void gsh_print(double, int *, FILE *);
void read_gsh_coeff(FILE *, int, int, int, double *, double *, double *, double *, double *, double *, int *);
/* interpolate_she_model.c */
int interpolate_she_model(struct mod *, struct mod *, double, int, int);
void mean_model(struct mod *, struct mod *, double, double, int, int);
double my_make_nan(void);
/* misc.c */
void mycalloc_cp(double **, int);
void mycalloc_dp(float **, int);
void mymalloc_dp(float **, int);
void mymalloc_cp(double **, int);
void myrealloc_dp(float **, int);
void myrealloc_cp(double **, int);
void zero_cp(double *, int);
void zero_dp(float *, int);
/* mod_modelbase.c */
void fit_base_functions(double, double, double *, double *, struct mod *, double (*)(double, double, double *, int, double), double (*)(int, int, int));
/* model2scatter.c */
/* modmodellmax.c */
/* mygrdio.c */
void grid_output(int, char *, float *, int, int, double, double, double, double, double, double, int, char **, int, unsigned short, void *);
void my_gmt_write_grd(float *, unsigned short, int, char **, char *, int, int, double, double, double, double, double, double, void *);
/* myopen.c */
FILE *myopen(char *, char *);
FILE *myopen_wn(char *, char *, char *);
/* nr_utils.c */
double gammln(double);
void nrerror(char []);
double *vector(long, long);
int *ivector(long, long);
unsigned char *cvector(long, long);
unsigned long *lvector(long, long);
double *dvector(long, long);
double **matrix(long, long, long, long);
double **dmatrix(long, long, long, long);
int **imatrix(long, long, long, long);
double **submatrix(double **, long, long, long, long, long, long);
double **convert_matrix(double *, long, long, long, long);
void free_vector(double *, long, long);
void free_ivector(int *, long, long);
void free_cvector(unsigned char *, long, long);
void free_lvector(unsigned long *, long, long);
void free_dvector(double *, long, long);
void free_matrix(double **, long, long, long, long);
void free_dmatrix(double **, long, long, long, long);
void free_imatrix(int **, long, long, long, long);
void free_submatrix(double **, long, long, long, long);
void free_convert_matrix(double **, long, long, long, long);
/* ones.c */
float checker(int, int, int, int, int);
/* plgndr.c */
void plgndr(double *, int, double *, int);
void plgndr_g(double *, int, int, double *);
void plgndr2(float *, int, double *, int);
double slgndr(int, int, double);
void pdtheta_lgndr(double *, int, double *, double *, int, double *, int);
/* plotlegendre.c */
/* powercorr.c */
double calc_correlation_model(struct mod *, int, int, int, int, int *);
double calc_total_power_model(struct mod *, int, int);
double calc_rms_model(struct mod *, int, int);
double degree_power_model(struct mod *, int, int);
double degree_power(double *, double *, int);
double degree_power_gsh(double *, double *, double *, double *, int);
double correlation_gsh(double *, double *, double *, double *, double *, double *, double *, double *, int, int, int, int, int);
void add_to_xy(double **, double **, int *, double, double);
double correlation(double *, double *, double *, double *, int, int, int, int);
double correlation_pt(double *, double *, double *, double *, double *, double *, double *, double *, int, int, int, int);
double ccl_correlation(double *, double *, int, int, int *, int);
double weighted_correlation(double *, double *, double *, double *, int, int, int);
double calc_rms(double *, double *, int);
double calc_rms_gsh(double *, double *, double *, double *, int);
double calc_total_power(double *, double *, int);
double calc_total_power_gsh(double *, double *, double *, double *, int);
/* rand.c */
double ran1(long *);
/* read_she_model.c */
void read_she_model(char *, struct mod *, int, int, int);
void allocate_model_coefficients(struct mod *);
void deallocate_model_coefficients(struct mod *);
void copy_model_par(struct mod *, struct mod *);
/* readflt.c */
/* scale_model.c */
/* select_lms.c */
void select_lms(int, int, int *, int *, int *);
/* shana.c */
void calc_coeff(float **, int, int, int, char *, int, float *, int, int, int, float, float, float, float, float, float, int, int, int, double);
void make_a(int, int, int, float **, float *, int, char *, int, int *, double, int, int, float *, int);
int aij(int, int, int, int, int);
double interpolate(double *, float *, double, int);
float dist_rad(float *, int, int);
void check_opmode(char *, int *, char *, int *, int *);
void gmt2myconvention_rotate4(float *, int, int, float);
void gmt2myconvention_rotate(float *, int, int, float, struct GMT_GRID *);
void phelp(char *);
void minmax(float *, float *, float *, float *, int);
void field_message(float, float, float, float, float, float, int, int, float, float, float, char *);
/* shsyn.c */
int nextpwrtwo(int);
void phelp(char *);
void read_ab(double *, double *, int, char *, FILE *, int, int, int);
void check_out_mode(char *, int *, char *);
/* spear.c */
double corr_sub(double *, double *, int, double *, int);
void pearsn(double [], double [], unsigned long, double *, double *, double *, int);
void spear(double [], double [], unsigned long, double *, double *, double *, double *, double *);
void crank(unsigned long, double [], double *);
void sort2(unsigned long, double [], double []);
double erfcc(double);
double betai(double, double, double);
double betacf(double, double, double);
/* sphex_lin_reg.c */
void sphex_lin_reg(double *, double *, double *, double *, double *, double *, int, int, int, double *);
/* spline.c */
void splint(double *, double *, double *, int, double, float *);
double *vector(long, long);
void free_vector(double *, long, long);
/* splinesc.c */
double spline_norm_damping(int, int, int);
/* test.c */
/* test_intel_stdio.c */
/* transpose.c */
/* write_coefficients.c */
void write_coefficients(double **, double **, int, double, int, FILE *, int);
void write_vector_coefficients(double **, double **, double **, double **, int, double, int, FILE *, int);
void write_model(struct mod *, double, int, FILE *);
void write_model_layer(struct mod *, int, double, int, int, FILE *);
