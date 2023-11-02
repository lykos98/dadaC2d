#include "./include/dadac.h"
#include <stdio.h>
#include <time.h>

struct Options 
{
	char* inputFile;
	char* outputFile;
    double Z;
    int halo;
    int k;
    int UseSparseBorders;
    int FileInFloat32;
    unsigned int data_dims;

};

void write_border_idx(const char * fname, Clusters * c)
{
    FILE * f = fopen(fname, "w");

    if(c -> UseSparseBorders)
    {
	    for(idx_t i = 0; i < c -> centers.count; ++i)
	    {
		for(idx_t el = 0; el < c ->SparseBorders[i].count; ++el)
		{
		    int a = c -> SparseBorders[i].data[el].idx; 
		    fprintf(f, "%d ",a);
		}
		fprintf(f,"\n");
	    }
    }
    else
    {
	    for(idx_t i = 0; i < c -> centers.count; ++i)
	    {
		for(idx_t j = 0; j < c -> centers.count; ++j)
		{
		    //int a = c -> borders[i][j].idx == NOBORDER ? -1 : c -> borders[i][j].idx;
		    //fprintf(f,"%d ",a);
		    int a = c -> borders[i][j].idx; 
		    if(c -> borders[i][j].idx != NOBORDER) fprintf(f, "%d ",a);
		}
		fprintf(f,"\n");
	    }
    }
    fclose(f);
}

void write_point_info(const char * fname, Datapoint_info * particles, idx_t n)
{
    FILE * f = fopen(fname,"w");
    for(idx_t i = 0; i < n; ++i)
    {
        fprintf(f,"%lu\t",(uint64_t)particles[i].kstar);
        fprintf(f,"%d\t", particles[i].cluster_idx);
	#ifdef USE_FLOAT32
		fprintf(f,"%.6f\t",particles[i].log_rho);
	#else
		fprintf(f,"%.11lf\t",particles[i].log_rho);
	#endif
        //fprintf(f,"%.12lf\t",particles[i].log_rho_c);
        //fprintf(f,"%.12lf\t",particles[i].log_rho_err);
        //fprintf(f,"%.12lf\t",particles[i].g);
        fprintf(f,"%d\t",particles[i].is_center);
        fprintf(f,"\n");
    }
    fclose(f);
}

void printHelp()
{
	printf("USAGE: ./driver i=[INPUT_FILE] o=[OUTPUT_FILE] d=[d] t=[t] z=[Z] h=[HALO] k=[k] s=[s] t=[t]\n");
	printf("\tINPUT_FILE : input file, file path\n");
	printf("\tOUTPUT_FILE: output file, file path\n");
	printf("\td	     : Lenght of the data vectors (int) \n");
	printf("\tZ	     : Z value, float\n");
	printf("\tHALO	     : Assign halo, y/n [yes/no]\n");
	printf("\tk	     : Number of neighbors to use, int (>0) \n");
	printf("\ts	     : Use sparse borders implementation, y/n [sparse/dense]\n");
	printf("\tt	     : Input binary is in Float32, y/n [float/double]\n");
	printf("\nThe program gives as output the cluster assignment of each datapoint\n");
	return;
}

void checkEq(char c)
{
	if(c != '='){
		printf("Wrongly formatted input\n");
		printHelp();
		exit(1);
	};
	return;
}

struct Options Parser(int argc, char** argv)
{
	struct Options opt = {.inputFile=NULL, .outputFile=NULL,.FileInFloat32 = 1,.k=1001, .halo = 1, .Z = 2, .data_dims = 0, .UseSparseBorders = 0 };
	if(argc < 2)
	{
		printHelp();
		exit(1);
	}
	for(int i = 1; i < argc; ++i)	
	{
		checkEq(argv[i][1]);
		switch(argv[i][0])
		{
			case 'i':
				opt.inputFile = argv[i] + 2;
				break;
			case 'o':
				opt.outputFile = argv[i] + 2;
				break;
			case 't':
				opt.FileInFloat32 = argv[i][2] == 'y' ? 1 : 0;
				break;
			case 'd':
				opt.data_dims = atoi(argv[i] + 2);
				break;
			case 'k':
				opt.k = atoi(argv[i] + 2);
				break;
			case 'h':
				opt.halo = argv[i][2] == 'y' ? 1 : 0;
				break;
			case 's':
				opt.UseSparseBorders = argv[i][2] == 'y' ? 1 : 0;
				break;
			case 'z':
				opt.Z = atof(argv[i] + 2);
				break;
			default:
				printHelp();
				exit(1);
				break;
		}
	}
	if(!(opt.inputFile) || !(opt.outputFile)){
		printf("Please provide input and output file paths\n");
		printHelp();
		exit(1);
	}
	if(opt.data_dims == 0)
	{
		printf("Please specify the lenght of each vector\n");
		printHelp();
		exit(1);
	}
	return opt;
}




int main(int argc, char** argv){

    //char aux_fname[80];

    //struct timespec start, finish;
    //double elapsed;
    struct timespec start_tot, finish_tot;
    double elapsed_tot;
    //Start timer
    clock_gettime(CLOCK_MONOTONIC, &start_tot);

	//struct Options opt = Parser(argc, argv);

    
    FILE* f = fopen("../euclid/img.dat","r");
    if(!f)
    {
        printf("Nope\n");
        exit(1);
    }
    fseek(f,0,SEEK_END);
    size_t n = ftell(f);
    rewind(f);

	int nrows = 7723;
	int ncols = 6945;
    n = n/(4);
    //n = nrows*ncols; 
    printf("Reading %lu particles\n",(uint64_t)n);


    FLOAT_TYPE* data = (FLOAT_TYPE*)malloc(n*sizeof(FLOAT_TYPE));

	if(1)
	{
		float* df = (float*)malloc(n*sizeof(float));
		size_t fff = fread(df,sizeof(float),n,f);
		//printf("Read %luB\n",fff);
		fclose(f);

		for(idx_t i = 0; i < n; ++i) data[i] = (FLOAT_TYPE)(df[i]);

		free(df);
	}
	else
	{
		double* df = (double*)malloc(n*sizeof(double));
		size_t fff = fread(df,sizeof(double),n,f);
		printf("Read %luB\n",fff);
		fclose(f);

		for(idx_t i = 0; i < n; ++i) data[i] = (FLOAT_TYPE)(df[i]);

		free(df);
	}

	//Datapoint_info* particles = NgbhSearch_kdtree(data, n, opt.data_dims, opt.k); 
	//Datapoint_info* particles = NgbhSearch_vptree(data, n,sizeof(FLOAT_TYPE), opt.data_dims, opt.k, eud); 
    /********************************
     * Intrinsic Dimension estimate *
     ********************************/

    //double id = idEstimate(particles, n);

    /***********************
     * Density computation *
     ***********************/
	int* mask = (int*)malloc(n*sizeof(int));
	FILE* ff = fopen("../euclid/mask.dat", "r");
	fread(mask, sizeof(int), n, ff);
	fclose(ff);
	//for(int i = 0; i < n; ++i) mask[i] = 1;
	float_t Z = 10;
	//Datapoint_info* computeDensityFromImg(FLOAT_TYPE* vals, int* mask, size_t nrows, size_t ncols);
    Datapoint_info* particles = computeDensityFromImg(data, mask, (int)nrows, (int)ncols, 15); 
    computeCorrection(particles,mask,n,Z);

    /********************
     * First clustering *
     ********************/

    Clusters c = Heuristic1(particles, mask, (int)nrows, (int)ncols);

    /***************************************************************************************
     * Allocate borders and other things to store clustering info                          *
     * Then Find borders between clusters and then merge clusters using peaks significance *
     ***************************************************************************************/
   // Clusters_allocate(&c);  
    Clusters_allocate(&c, 1);  

    // sprintf(aux_fname, "%s_int", argv[2]);
    // write_point_info(aux_fname,particles,n);

    Heuristic2(&c, particles, mask, (int)nrows, (int)ncols);

    //sprintf(aux_fname, "%s_bord_int", argv[2]);
    //write_border_idx(aux_fname,&c);

    c.n = n;
    
    Heuristic3(&c, particles, Z, 1);

	write_point_info("../euclid/out.dat", particles, n);
	free(data);
	freeDatapointArray(particles,n);
	Clusters_free(&c);


    clock_gettime(CLOCK_MONOTONIC, &finish_tot);
    elapsed_tot = (finish_tot.tv_sec - start_tot.tv_sec);
    elapsed_tot += (finish_tot.tv_nsec - start_tot.tv_nsec) / 1000000000.0;
    printf("ELAPSED time (measured by driver): %.3lfs\n\n", elapsed_tot);

    return 0;
}
