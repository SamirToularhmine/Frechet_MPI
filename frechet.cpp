#include <iostream>
#include <iostream>
#include <fstream>
#include <string>
#include <math.h>
#include <vector>
#include <mpi.h>

using namespace std;

double* lecture (const char* filename, int* n)
{
    ifstream f(filename);
    int size;
    f>>size;
    (*n) = size;
    double* t = new double[2*size];
    string line;
    int ptr = 0;
    f>>line;
    while (!f.eof()) {
        int pos = line.find(";");
        t[ptr] = stod(line.substr(0, pos));
        ptr++;
        t[ptr] = stod(line.substr(pos + 1));
        ptr++;
        f>>line;
    }
    f.close();
    return t;
}

double dist_e(double* v1, double* v2) {
    return sqrt((v2[0]-v1[0])*(v2[0]-v1[0]) + (v2[1]-v1[1])*(v2[1]-v1[1]));
}

double* matrice_distance(double* traj1, int n1, double* traj2, int n2) {

    /** Matrice de distances optimisée **/
    double max_diag = 0;
    double max_diagl_global;
    int min = std::min(n1, n2);
    int max = std::max(n1, n2);

    int pid, nprocs;  
    MPI_Status status;
    MPI_Comm_rank(MPI_COMM_WORLD, &pid);
    MPI_Comm_size (MPI_COMM_WORLD, &nprocs);

    // Ici on trouve la position actuelle dans le tableau traj2, sinon on ne peux pas commencer au
    // bon endroit dans traj1 pour calculer la diagonale. Merci Ionas.

    int reste = n2 % nprocs;
    int ilocal = pid < reste ? (n2 / nprocs) * pid + pid : (n2 / nprocs) * pid + reste;
    int tailleL2 = n2 / nprocs;
    if(pid < reste){
        tailleL2++;
    }

    // On initialise distance en fonction du nombre d'element dans le morceau de L2 reçu.
    double* distances = new double[n1 * tailleL2];

    // Calcul de la diagonale + de son max.
    for (int i=0; i < tailleL2; i++){
        double curr_dist = dist_e(traj1+(2 * (ilocal + i)), traj2+(2*i));
        distances[(ilocal + i + (i * n1))] = curr_dist;

        if(curr_dist > max_diag){
            max_diag = curr_dist;
        }
    }

    // Si la trajectoire L1 a + d'éléments que la trajectoire L2, alors il faut calculer la fausse diagonale restante.
    // On suppose donc que si le nombre d'élement des deux trajectoires est différents, alors L1 devra avoir le plus grand nombre.
    if(pid == (nprocs - 1) && n1 > n2){
        for(unsigned int i = 0; i < n1 - n2; i++){
            double curr_dist = dist_e(traj1+(2 * (n2 + i)), traj2+(2*(tailleL2 - 1)));
            distances[(ilocal + (i + 1))] = curr_dist;

            if(curr_dist > max_diag){
                max_diag = curr_dist;
            }

        }
    }

    // On applique une opération de réduction pour trouver le max sur la diagonale.
    MPI_Allreduce(&max_diag, &max_diagl_global, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

    // Parcours triangulaire inférieur + calcul de la valeur + insertion dans distances.
    for(int i = 0; i < tailleL2; i++){
        bool colBigger = false;
        for(int j = ilocal + 1; j < n1; j++){
            if(colBigger){
                distances[(i * n1) + j] = INFINITY;
            }else{
                double curr_dist = dist_e(traj1+(2*j), traj2+(2*i));
                if(curr_dist > max_diagl_global){
                    distances[(i * n1) + j] = INFINITY;
                    colBigger = true;
                }else{
                    distances[(i * n1) + j] = curr_dist;
                }
            }
        }
    }

    // Parcours triangulaire supérieur + calcul de la valeur + insertion dans distances.
    for(int i = 0; i < ilocal + tailleL2 - 1; i++){
        bool rowBigger = false;
        for(int j = 0; j < tailleL2; j++){
            if(i == ilocal + j){
                continue;
            }
            if(rowBigger){
                distances[(j * n1) + i] = INFINITY;  
            }else{
                double curr_dist = dist_e(traj1+(2*i), traj2+(2*j));
                if(curr_dist > max_diagl_global){
                    rowBigger = true;
                    distances[(j * n1) + i] = INFINITY;
                }else{
                    distances[(j * n1) + i] = curr_dist;
                }
            }
            
        }
    }

    return distances;
}


double min(double a, double b) {
    return ((a<b)?a:b);
}
double max(double a, double b) {
    return ((a<b)?b:a);
}

double* matrice_frechet(double* distances, int n1, int n2)
{
    double* frechet = new double[n1*n2];

    frechet[0] = distances[0];

    for (int i=1; i<n1; i++)
        frechet[i*n2] = max(distances[i*n2],distances[(i-1)*n2]);

    for (int j=1; j<n2; j++)
        frechet[j] = max(distances[j], distances[(j-1)*n2]);

    for (int i=1; i<n1; i++)
        for (int j=1; j<n2; j++)
            frechet[i*n2+j] = max(distances[i*n2+j],min(frechet[(i-1)*n2+j],min(frechet[i*n2+(j-1)],frechet[(i-1)*n2+(j-1)])));

    return frechet;

}

int main(int argc, char*argv[]) {
    int* n_traj = new int[2];
    double* traj1;
    double* traj2;
    int* sendcounts;
    int* displs;
    double* traj2_scattered;

    int pid, nprocs;  
    MPI_Status status;
    MPI_Init (&argc , &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &pid);
    MPI_Comm_size (MPI_COMM_WORLD, &nprocs);

    if(pid == 0){
        // On ne lit que les trajectoires dans le processeur 0
        traj1 = lecture(argv[1],&n_traj[0]);
        traj2 = lecture(argv[2],&n_traj[1]);

        // On initialise displs et sendcounts uniquement dans le processeur 0
        sendcounts = new int[nprocs];
        displs = new int[nprocs];
        int reste = n_traj[1] % nprocs;

        for(unsigned int i = 0; i < nprocs; i++){
            sendcounts[i] = 2 * (n_traj[1] / nprocs);
            if(reste > 0){
                sendcounts[i] += 2;
                reste --;
            }
        }

        displs[0] = 0;
        for(int i = 1; i < nprocs; i++){
            displs[i] = sendcounts[i-1] + displs[i-1];
        }
    
        for (int i=0; i<n_traj[0]; i++)
            cout << traj1[2*i] << " " << traj1[2*i+1] << endl;
        for (int i=0; i<n_traj[1]; i++)
            cout << traj2[2*i] << " " << traj2[2*i+1] << endl;

    }

    // On broadcast le nombre d'élements des tableaux de distances.
    MPI_Bcast(n_traj, 2, MPI_INT, 0, MPI_COMM_WORLD);

    if(pid != 0){
        traj1 = new double[2 * n_traj[0]];
    }

    // On broadcast la trajectoire L1.
    MPI_Bcast(traj1, 2 * n_traj[0], MPI_DOUBLE, 0, MPI_COMM_WORLD);

    int nlocal = n_traj[1] / nprocs;
    int reste = n_traj[1] % nprocs;

    if(pid < reste){
        nlocal ++;
    }

    nlocal *= 2;

    // traj2_scattered correspond au morceau de L2 reçu pour chaque processeur.
    traj2_scattered = new double[nlocal];

    // On découpe L2 en fonction du nombre de processeurs.
    MPI_Scatterv(traj2, sendcounts, displs, MPI_DOUBLE, traj2_scattered, nlocal, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    // On calcule la matrice de distance pour chaque processeur.
    double* mat_dist = matrice_distance(traj1,n_traj[0],traj2_scattered,n_traj[1]);

    double* distances;
    int* distances_displs;
    int* distances_recvcounts;
    if(pid == 0){
        distances = new double[n_traj[0] * n_traj[1]];
        distances_displs = new int[nprocs];
        distances_recvcounts = new int[nprocs];

        distances_displs[0] = 0;
        for(unsigned int i = 1; i < nprocs; i++){
            distances_displs[i] = (n_traj[0] * (sendcounts[i-1] / 2)) + distances_displs[i-1];
        }
        for(int i = 0; i < nprocs; i++){
            distances_recvcounts[i] = (sendcounts[i] / 2) * n_traj[0];
        }
    }

    // On rassemble les tableaux de distances des différents processeurs dans le tableau "distances".
    MPI_Gatherv(mat_dist, (n_traj[0] * (nlocal / 2)) , MPI_DOUBLE, distances, distances_recvcounts, distances_displs, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    double* distances_transposee;

    // Comme on a rassembler à la suite les différentes colonnes de la matrice de distance finale,
    // Il faut transposer ce tableau afin que ces colonnes deviennent vraiment des colonnes et non pas
    // des lignes comme avant la transposition. Merci Mme.Robert.
    if(pid == 0){
        // Transposition du tableau
        distances_transposee = new double[n_traj[0] * n_traj[1]];
        for(int i = 0; i < n_traj[0]; i++){
            for(int j = 0; j < n_traj[1]; j++){
                distances_transposee[i * n_traj[1] + j] = distances[i + (j * n_traj[0])];
            }
        }
    }

    if(pid == 0){
        cout << "la matrice de distance" << endl;

        for (int i=0; i<n_traj[0]; i++) {
            for (int j = 0; j < n_traj[1]; j++)
                cout << distances_transposee[i * n_traj[1] + j] << "  ";
            cout << endl;
        }

        double* frechet = matrice_frechet(distances_transposee,n_traj[0],n_traj[1]);

        cout << "la matrice de fréchet" << endl;

        for (int i=0; i<n_traj[0]; i++) {
            for (int j = 0; j < n_traj[1]; j++)
                cout << frechet[i * n_traj[1] + j] << " ";
            cout << endl;
        }
    }
    
    MPI_Finalize();

    // Libérations de mémoire.
    delete[] traj2_scattered;

    if(pid == 0){
        delete[] displs;
        delete[] sendcounts;
        delete[] distances;
        delete[] distances_displs;
        delete[] distances_transposee;
    }

    delete[] mat_dist;

    return 0;
}
