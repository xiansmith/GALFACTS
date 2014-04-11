

void compute_hanning_coefficients(double H[], int n);

void compute_boxcar_coefficients(double H[], int n);

void interpolate_missing_data(double A[], int size);

void apply_smoothing_function(const double A[], double B[], int size, const double H[], int n);

void hanning_smooth_data(double A[], int size, int n);
