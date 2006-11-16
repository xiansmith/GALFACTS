
determine_correction_factors: determine_correction_factors.c
	gcc  -g -lm determine_correction_factors.c -o determine_correction_factors

leakage_examine: leakage_examine.c fetch_from_best_fit_model.o
	gcc  -g -lg2c -lm fetch_from_best_fit_model.o leakage_examine.c -o leakage_examine 

leakagetest: leakagetest.c fetch_from_best_fit_model.o
	gcc  -g -lg2c -lm fetch_from_best_fit_model.o leakagetest.c -o leakagetest
	
apply_correction_factors: apply_correction_factors.c fetch_from_best_fit_model.o
	gcc  -g -lg2c -lm fetch_from_best_fit_model.o apply_correction_factors.c -o apply_correction_factors

fetch_from_best_fit_model.o: fetch_from_best_fit_model.f
	g77 -g -c fetch_from_best_fit_model.f

clean:
	rm -f *.o leakage
