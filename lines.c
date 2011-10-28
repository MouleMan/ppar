#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <mpi.h>
#include <omp.h>
#include <errno.h>


#define EPSILON 1.0e-5
#define REAL_T float
#define C 0.2F
#define STEP 40

#define TAG_MSG_UP 0
#define TAG_MSG_DOWN 1

#define TAG_TO_STRING(TAG) ((TAG) == 0)?"TAG_MSG_UP":"TAG_MSG_DOWN"

#define error(MESSAGE, USE_ERROR) (err(__FILE__, __LINE__, (MESSAGE), (USE_ERROR)))

typedef struct proc
{
    int size;
    int rank;
    int gridW, gridH;
    int w, h;
    REAL_T **buff[2];
    REAL_T *buffHaut;
    REAL_T *buffBas;
    int iter;
    REAL_T maxError;
} proc;


/* Prototypes */
double my_gettimeofday();
void printGrid(FILE *f, proc *my, int id);

void err(const char * file, int line, char * message, int usePerror)
{
	if(usePerror)
		fprintf(stderr, "Fichier : \"%s\", ligne \"%d\" Message : \"%s\" Erreur : \"%s\" \n", 
				file, line, message, strerror(errno));
	else
		fprintf(stderr, "Fichier : \"%s\", à ligne \"%d\" Message : \"%s\"\n", 
				file, line, message);
	
    MPI_Finalize();
	exit(EXIT_FAILURE);
}

/**
 * Initalise le processus à partir des arguments d'entrée.
 * 
 * \return -1 si la fonction échoue, 0 sinon.
 */
static int initProc(proc * my, int * argc, char *** argv)
{
    int i, j;
    int start;

    memset(my, 0, sizeof(proc));

    MPI_Init(argc, argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &my->rank);
    MPI_Comm_size(MPI_COMM_WORLD, &my->size);

    my->gridW = 1000;
    my->gridH = 1000;

    if(*argc > 1)
        my->gridW = atoi((*argv)[1]);
    if(*argc > 2)
        my->gridH = atoi((*argv)[2]);

    if(my->rank == 0)
        printf("Dims: (%d * %d) Nb points: %i Taille memoire: %fMo X2\n", 
                my->gridW, my->gridH, my->gridW * my->gridH, (float) my->gridW * my->gridH * sizeof(REAL_T) / (1024*1024));

    if(my->gridH % my->size)
    {
        fprintf(stderr, "Le nombre de processus n'est pas un multiple de la hauteur de l'image\n");
        return EXIT_FAILURE;
    }

    if(my->rank != 0)
    {
        my->buffHaut = malloc(sizeof(REAL_T*) * my->gridW);
        if(my->buffHaut == NULL)
            error("Erreur d'allocation du buffer", 1);
    }
    
    if(my->rank != my->size - 1)
    {
        my->buffBas = malloc(sizeof(REAL_T*) * my->gridW);
        if(my->buffBas == NULL)
            error("Erreur d'allocation du buffer", 1);
    }

    my->w = my->gridW;
    my->h = my->gridH / my->size;
                
    for(i = 0; i < 2; i++)
    {
        my->buff[i] = malloc(sizeof(REAL_T*) * my->gridH);

        if(my->buff[i] == NULL)
            error("Erreur d'allocation du buffer", 1);

        for(j = 0; j < my->gridH; j++)
        {
            my->buff[i][j] = malloc(sizeof(REAL_T) * my->gridW);
            if(my->buff[i][j] == NULL)
                error("Erreur d'allocation du buffer", 1);
        }
    }

    for(i = 0; i < my->h; i++)
        for(j = 0; j < my->w; j++)
            my->buff[0][i][j] = my->rank;

    /* Start est le numero de la ligne, dans le buffer local,
     * ou se trouve le premier STEP.
     * Il faut prendre en compte le rang, ainsi que la hauteur
     * variable des buffers
     */
     /*
    start = (int)ceil(my->h * my->rank / (float)STEP);
    start = (start * STEP) - my->rank * my->h;
    
    for(i = start; i < my->gridH; i += STEP)
    {
 */       /*printf("%d i = %d\n", my->rank, i);*/
   /*     for(j = 0; j < my->w; j++)
            my->buff[0][i][j] = (REAL_T) 100.0;
    }*/
    
    
    
    
//    printf("proc %d: (%d * %d) Nb points: %i Taille memoire: %fMo X2\n", 
                //my->rank, my->w, my->h, my->w * my->h, (float) my->w * my->h * sizeof(REAL_T) / (1024*1024));
//
    //printGrid(my->rank ? stdout : stderr, my, 0);
    
    return 0;
}


int compute(proc * my)
{
    int i, j;
    
    MPI_Request recvRequest[2];
    MPI_Request sendRequest[2];
    MPI_Status status[2];
    int cptRequest = 0;
    
    do
    {
        REAL_T **prevBuff = my->buff[my->iter % 2];
        REAL_T **currBuff = my->buff[(my->iter + 1) % 2];

        /* EMISSION DES BUFFERS DE L'ITERATION PRECEDENTE */
        
        /* Processus différent de 0 */
        if(my->rank > 0)
        {
            // multiplier les Isend() pour la même ligne, pour commencer à calculer avant de tout avoir 
            
            /* Emission de la ligne du haut au processus du dessus */
            MPI_Isend(prevBuff[0], my->w, MPI_FLOAT, my->rank - 1, TAG_MSG_UP, MPI_COMM_WORLD, &sendRequest[cptRequest]);
            
            printf("%d envoie une requête TAG_MSG_UP a %d \n", my->rank, my->rank - 1);
            printf("%d attend une requête TAG_MSG_DOWN de %d \n", my->rank, my->rank - 1);
            
            /* Réception la ligne du haut depuis le processus du dessus */
            MPI_Irecv(my->buffHaut, my->w, MPI_FLOAT, my->rank - 1, TAG_MSG_DOWN, MPI_COMM_WORLD, &recvRequest[cptRequest]);
            
            cptRequest++;
        }
        
        /* Processus différent de N */
        if(my->rank < my->size - 1)
        {
            /* Emission de la ligne du bas au processus du dessous */
            MPI_Isend(prevBuff[my->h - 1], my->w, MPI_FLOAT, my->rank + 1, TAG_MSG_DOWN, MPI_COMM_WORLD, &sendRequest[cptRequest]);
            
            printf("%d envoie une requête TAG_MSG_DOWN a %d \n", my->rank, my->rank + 1);
            printf("%d attend une requête TAG_MSG_UP de %d \n", my->rank, my->rank + 1);
            
            /* Réception la ligne du bas depuis le processus du dessous */
            MPI_Irecv(my->buffBas, my->w, MPI_FLOAT, my->rank + 1, TAG_MSG_UP, MPI_COMM_WORLD, &recvRequest[cptRequest]);
            
            cptRequest++;
        }
        
        /* Vérification des réceptions */
        for(i = 0; i < cptRequest; i++)
        {
            /* TODO dangereux la boucle ? */
            int index;
            MPI_Waitany(cptRequest, recvRequest, &index, &status[i]);    
            printf("%d recptionne une requête \"%s\" de \"%d\" 1ere case %f\n", 
                my->rank, TAG_TO_STRING(status[i].MPI_TAG), status[i].MPI_SOURCE, (status[i].MPI_SOURCE == my->rank - 1) ? *my->buffHaut : *my->buffBas);
            
        }
        
        /* Vérifications des envois */
        MPI_Waitall(cptRequest, sendRequest, status);
        printf("%d Tout a été envoyé \n", my->rank);
        
        while(1);
        
        
        
        
        /*
        for(i = 0; i < my->gridH; i++)
        {
            for(j = 0; j < my->gridW; j++)
            {
                currBuff[i][j] = (1.0F - 4.0F * C) * prevBuff[i][j]
                    + C * (j > 0   ? prevBuff[i][j - 1] : prevBuff[i][j])
                    + C * (j < my->gridW-1 ? prevBuff[i][j + 1] : prevBuff[i][j])
                    + C * (i > 0   ? prevBuff[i - 1][j] : prevBuff[i][j]) 
                    + C * (i < my->gridH-1 ? prevBuff[i + 1][j] : prevBuff[i][j]);
            }
        }*/

        /* calcul de l'erreur */
        /*
        my->maxError = 0.0;

        // /!\ TODO: on ne doit parcourir que les points calculés par nous à l'itération courrante
        for(i = 0; i < my->h; i++)
        {
            for(j = 0; j < my->w; j++)
            {
                REAL_T error = (REAL_T) fabs(currBuff[i][j] - prevBuff[i][j]);

                if(error > my->maxError)
                    my->maxError = error; 
            }
        }
        */
        //fprintf(stderr, "Iteration %d : \t delta = %.3e\n", iter, maxError);

        my->iter++;
    }
    while(my->maxError > EPSILON);

    return (my->iter + 1) % 2;
}


void printGrid(FILE *f, proc *my, int id)
{
    int i, j;

    fprintf(f, "===%d==\n", my->rank);

    for(i = 0; i < my->h; i++)
    {
        for(j = 0; j < my->w; j++)
            fprintf(f, "%.2f ", my->buff[id][i][j]);
            
        fprintf(f, "\n");
    }
}


int main(int argc, char **argv)
{
    double debut, fin;
    proc my;

    initProc(&my, &argc, &argv);
    
    /* Attente de tout le monde pour lancer le calcul */
    MPI_Barrier(MPI_COMM_WORLD);

    if(my.rank == 0)
        debut = my_gettimeofday();

    compute(&my);

    /* Attente de tout le monde pour arretter le chrono */
    MPI_Barrier(MPI_COMM_WORLD);

    if(my.rank == 0)
    {
        /* Arret horloge : */
        fin = my_gettimeofday();

        /* Affichage des resultats : */
        printf("Nombre d'iterations : %d\n", my.iter);
        printf("Erreur (delta) = %.3e\n", my.maxError);
        printf("Temps total : %.1f s\n", fin - debut);
    }

    free(my.buff[0]);
    free(my.buff[1]);

    MPI_Finalize();

    return EXIT_SUCCESS;
}


/* Gestion du temps Linux / Visual Studio */
#ifdef _MSC_VER

#include <sys/timeb.h>

double my_gettimeofday()
{
    struct _timeb timebuffer;
    _ftime (&timebuffer);
    return timebuffer.time + timebuffer.millitm * 1.0e-3;
}


#else

#include <sys/time.h>

double my_gettimeofday()
{
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return tv.tv_sec + tv.tv_usec * 1.0e-6;
}

#endif
