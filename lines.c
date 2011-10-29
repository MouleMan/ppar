#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <mpi.h>
#include <omp.h>
#include <errno.h>


#define EPSILON 1.0e-3
#define REAL_T float
#define C 0.2F
#define STEP 100

#define TAG_MSG_UP 0
#define TAG_MSG_DOWN 1

#define TAG_TO_STRING(TAG) ((TAG) == 0)?"TAG_MSG_UP":"TAG_MSG_DOWN"

#define error(MESSAGE, USE_ERROR) (err(__FILE__, __LINE__, (MESSAGE), (USE_ERROR)))



typedef struct proc
{
    /* Données générales sur l'execution*/
    int size;                       // nombre de processus
    int rank;                       // rang du processus
    int gridW, gridH;               // taille totale de la grille

    /* Sorties du programme */
    int iter;                       // nombre d'itérations
    REAL_T delta;                   // delta de la dernière itération

    /* Données sur les buffers locaux au processus */
    int w, h;                       // taille des buffers locaux
    REAL_T **buff[2];               // buffers locaux
    REAL_T *buffHaut;               // ligne calculée par le processus du dessus
    REAL_T *buffBas;                // ligne calculée par le processus du dessous
    REAL_T **prevBuff;              // pointeur permettant de swapper les buffers d'une itération sur l'autre
    REAL_T **currBuff;              // idem

    /* Données annexes */
    MPI_Request recvRequest[2];     // requêtes MPI sur les réceptions
    MPI_Request sendRequest[2];     // requêtes MPI sur les émissions
    int cptRequest;                 // nombre de requêtes MPI effectuées
} proc;


/* Prototypes */
double my_gettimeofday();
void printGrid(FILE *f, proc *my, int id);
void err(const char * file, int line, char * message, int usePerror);
static int initProc(proc * my, int * argc, char *** argv);                                                                                           // static => \TODO Trouver un autre binôme qui ne pinaille pas au point d'essayer d'empecher qu'une fonction puisse être appellée depuis un autre fichier alors que le projet n'en comporte qu'un seul et que, par conséquent, c'est pas près d'arriver
void freeProc(proc *my);
void run(proc * my);
void swapUp(proc * my);
void swapDown(proc * my);
void computeCenter(proc *my);
void computeBorders(proc *my);
REAL_T computeDelta(proc *my);


int main(int argc, char **argv)
{
    double debut, fin;
    proc my;

    initProc(&my, &argc, &argv);

    /* Attente de tout le monde pour lancer le calcul */
    MPI_Barrier(MPI_COMM_WORLD);
    debut = my_gettimeofday();

    /* Calcul parallélisé */
    run(&my);

    /* Attente de tout le monde pour arretter le chrono */
    MPI_Barrier(MPI_COMM_WORLD);
    fin = my_gettimeofday();

    /* Affichage des resultats */
    if(my.rank == 0)
    {
        printf("Nombre d'iterations : %d\n", my.iter);
        printf("Erreur (delta) = %.3e\n", my.delta);
        printf("Temps total : %.1f s\n", fin - debut);
    }

    freeProc(&my);

    MPI_Finalize();

    /* Empeche la console de se fermer à la fin de l'execution sous windows */
#ifdef _WIN32
    getchar();
#endif

    return EXIT_SUCCESS;
}


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

    /* Affichage d'informations globales */
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


    /* Initialisation du buffer */
    for(i = 0; i < my->h; i++)
        for(j = 0; j < my->w; j++)
            my->buff[0][i][j] = (REAL_T) 25.0;

    /* Start est le numero de la ligne, dans le buffer local,
    * ou se trouve le premier STEP.
    * Il faut prendre en compte le rang, ainsi que la hauteur
    * variable des buffers */
    start = (int)ceil(my->h * my->rank / (float)STEP);
    start = (start * STEP) - my->rank * my->h;

    for(i = start; i < my->gridH; i += STEP)
    {
        //printf("%d i = %d\n", my->rank, i);
        for(j = 0; j < my->w; j++)
            my->buff[0][i][j] = (REAL_T) 100.0;
    }
    
    /*
    for(i = 0; i < my->h; i++)
        for(j = 0; j < my->w; j++)
            my->buff[0][i][j] = my->rank * 10000.F + i;
    */

    /* Affichage d'informations locales au processus */
    printf("proc %d: (%d * %d) Nb points: %i Taille memoire: %fMo X2\n", 
        my->rank, my->w, my->h, my->w * my->h, (float) my->w * my->h * sizeof(REAL_T) / (1024*1024));

    //printGrid(my->rank ? stdout : stderr, my, 0);

    return 0;
}


void freeProc(proc *my)
{
    int i, j;

    if(my->buffBas)
        free(my->buffBas);

    if(my->buffHaut)
        free(my->buffHaut);

    for(i = 0; i < 2; i++)
    {
        for(j = 0; j < my->gridH; j++)
            free(my->buff[i][j]);

        free(my->buff[i]);
    }
}


void run(proc * my)
{
    MPI_Status status[2];
    REAL_T delta;

    do
    {
        /* Swap des buffers */
        my->prevBuff = my->buff[my->iter % 2];
        my->currBuff = my->buff[(my->iter + 1) % 2];
        my->cptRequest = 0;

        /* Demarage des echanges de données non bloquants avec les voisins */
        swapUp(my);
        swapDown(my);

        /* Calcul des valeurs de la grille pour la nouvelle itération */
        computeCenter(my);

        /* Idem pour les lignes dépendantes des voisins */
        computeBorders(my);

        /* Calcul du delta local de l'itération */
        delta = computeDelta(my);

        /* Attente de finalisation des emissions */
        MPI_Waitall(my->cptRequest, my->sendRequest, status);
        MPI_Allreduce(&delta, &my->delta, 1, MPI_FLOAT, MPI_MAX, MPI_COMM_WORLD); // TODO remplace le WaitAll() je pense

        //fprintf(stderr, "%d === Fin de l'iteration %d Delta = %.3e ====\n", my->rank, my->iter, my->delta);

        my->iter++;
    }
    while(my->delta > EPSILON);
}


/* Echange de lignes avec le processus du dessus */
void swapUp(proc * my)
{
    if(my->rank > 0)
    {
        //TODO multiplier les Isend() pour la même ligne, pour commencer à calculer avant de tout avoir 

        /* Emission de la ligne du haut au processus du dessus */
        MPI_Isend(my->prevBuff[0], my->w, MPI_FLOAT, my->rank - 1, TAG_MSG_UP, MPI_COMM_WORLD, &my->sendRequest[my->cptRequest]);

        //printf("%d envoie une requete TAG_MSG_UP a %d \n", my->rank, my->rank - 1);
        //printf("%d attend une requete TAG_MSG_DOWN de %d \n", my->rank, my->rank - 1);

        /* Réception la ligne du haut depuis le processus du dessus */
        MPI_Irecv(my->buffHaut, my->w, MPI_FLOAT, my->rank - 1, TAG_MSG_DOWN, MPI_COMM_WORLD, &my->recvRequest[my->cptRequest]);

        my->cptRequest++;
    }
}


/* Echange de lignes avec le processus du dessous */
void swapDown(proc * my)
{
    if(my->rank < my->size - 1)
    {
        /* Emission de la ligne du bas au processus du dessous */
        MPI_Isend(my->prevBuff[my->h - 1], my->w, MPI_FLOAT, my->rank + 1, TAG_MSG_DOWN, MPI_COMM_WORLD, &my->sendRequest[my->cptRequest]);

        //printf("%d envoie une requete TAG_MSG_DOWN a %d \n", my->rank, my->rank + 1);
        ///printf("%d attend une requete TAG_MSG_UP de %d \n", my->rank, my->rank + 1);

        /* Réception la ligne du bas depuis le processus du dessous */
        MPI_Irecv(my->buffBas, my->w, MPI_FLOAT, my->rank + 1, TAG_MSG_UP, MPI_COMM_WORLD, &my->recvRequest[my->cptRequest]);

        my->cptRequest++;
    }
}


REAL_T computeDelta(proc *my)
{
    int i, j;
    REAL_T delta = 0;

    for(i = 0; i < my->h; i++)
    {
        for(j = 0; j < my->w; j++)
        {
            REAL_T error = (REAL_T) fabs(my->currBuff[i][j] - my->prevBuff[i][j]);

            if(error > delta)
                delta = error; 
        }
    }

    return delta;
}


void computeBorders(proc *my)
{
    REAL_T *lineAbove, *lineUnder;
    MPI_Status status[2];
    int index;
    int i, j;

    /* Vérification des réceptions */
    for(i = 0; i < my->cptRequest; i++)
    {
        /* TODO dangereux la boucle ? */

        MPI_Waitany(my->cptRequest, my->recvRequest, &index, &status[i]);
        //printf("%d recptionne une requete \"%s\" de \"%d\" 1ere case %.2f\n", 
        //    my->rank, TAG_TO_STRING(status[i].MPI_TAG), status[i].MPI_SOURCE, (status[i].MPI_TAG == TAG_MSG_UP) ? my->buffBas[0] : my->buffHaut[0]);

        /* Le numero de La ligne qu'on est maintenant en mesure de calculer */
        index = (status[i].MPI_TAG == TAG_MSG_UP) ? my->h - 1 : 0;
        /* Le pointeur sur la ligne du dessus */
        lineAbove = (status[i].MPI_TAG == TAG_MSG_UP) ? my->prevBuff[index - 1] : my->buffHaut;
        /* Le pointeur sur la ligne du dessous */
        lineUnder = (status[i].MPI_TAG == TAG_MSG_UP) ? my->buffBas : my->prevBuff[index + 1];

        //printf("%d calcule sa ligne %d avec les valeurs haut %.2f et bas %.2f\n", my->rank, index, lineAbove[0], lineUnder[0]);

        /* Calcul de l'itération suivante pour la ligne reçue */
        for(j = 0; j < my->gridW; j++)
            my->currBuff[index][j] = (1.0F - 4.0F * C) * my->prevBuff[index][j]
            + C * (j > 0   ? my->prevBuff[index][j - 1] : my->prevBuff[index][j])
            + C * (j < my->gridW-1 ? my->prevBuff[index][j + 1] : my->prevBuff[index][j])
            + C * lineAbove[j] 
            + C * lineUnder[j];
    }
}


void computeCenter(proc *my)
{
    int i, j;

    /* Calcul de l'itération suivante pour tout sauf les lignes adjacentes à 
    * un voisin (premiere et/ou derniere selon le rang) */
    for(i = (my->rank != 0); i < my->gridH - (my->rank != my->size - 1); i++)
    {
        // TODO: traiter les premières et dernières colones en dehors de la boucle pour voir les perfs... Utilisable dans la boucle plus bas également??
        for(j = 0; j < my->gridW; j++)
        {
            my->currBuff[i][j] = (1.0F - 4.0F * C) * my->prevBuff[i][j]
                + C * (j > 0   ? my->prevBuff[i][j - 1] : my->prevBuff[i][j])
                + C * (j < my->gridW-1 ? my->prevBuff[i][j + 1] : my->prevBuff[i][j])
                + C * (i > 0   ? my->prevBuff[i - 1][j] : my->prevBuff[i][j]) 
                + C * (i < my->gridH-1 ? my->prevBuff[i + 1][j] : my->prevBuff[i][j]);
        }
    }
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


/* Gestion du temps Visual Studio */
#ifdef _MSC_VER

#include <sys/timeb.h>

double my_gettimeofday()
{
    struct _timeb timebuffer;
    _ftime (&timebuffer);
    return timebuffer.time + timebuffer.millitm * 1.0e-3;
}


/* Gestion du temps linux */
#else

#include <sys/time.h>

double my_gettimeofday()
{
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return tv.tv_sec + tv.tv_usec * 1.0e-6;
}

#endif
