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

    double* distances = new double[n1 * tailleL2];

    for(unsigned int i = 0; i < n1 * tailleL2; i++){
        distances[i] = 0;
    }

    // Calcul de la diagonale + de son max.
    for (int i=0; i < tailleL2; i++){
        double curr_dist = dist_e(traj1+(2 * (ilocal + i)), traj2+(2*i));
        distances[(ilocal + i + (i * n1))] = curr_dist;
        if(curr_dist > max_diag){
            max_diag = curr_dist;
        }
    }

    if(pid == (nprocs - 1) && n1 > n2){
        for(unsigned int i = 0; i < n1 - n2; i++){
            double curr_dist = dist_e(traj1+(2 * i), traj2+(2*(tailleL2 - 1)));
            distances[(ilocal + (i + 1))] = curr_dist;
            if(curr_dist > max_diag){
                max_diag = curr_dist;
            }
        }
    }
    //std::cout << " pid : " << pid << std::endl;

    /*for(int i = 0; i < n1 * tailleL2; i++){
        std::cout << distances[i] << " ";
    }
    std::cout << std::endl;*/

    MPI_Allreduce(&max_diag, &max_diagl_global, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

    /*int count = 0;
    double max_local = 0;
    std::vector<double> diags;
    for(int i = 0; i < max; i+= nprocs){
        if(i + pid < max){
            double curr_dist = 0;
            if(i + pid >= min){
                if(n1 < n2){
                    curr_dist = dist_e(traj1+2*(min - 1), traj2+2*(i + pid));
                }else{
                    curr_dist = dist_e(traj1+2*(i + pid), traj2+2*(min - 1));
                }
            }else{
                curr_dist = dist_e(traj1+2*(i + pid), traj2+2*(i + pid));
            }
            if(curr_dist > max_local){
                max_local = curr_dist;
            }
            diags.push_back(curr_dist);
            count++;
        }
    }*/

    // J'ai la diagonale + le max de la diagonale

    if(pid == 0){
        for(int i = 0; i < max; i++){
            //distances[(pid * i * min) + (i * nprocs)] = buff[i];
        }
    }

    // Calcul de la diagonale + de son max.
    /* for (int i=0; i<min; i++){
        double curr_dist = dist_e(traj1+2*i, traj2+2*i);
        distances[i*min+i] = curr_dist;
        if(curr_dist > max_diag){
            max_diag = curr_dist;
        }
    }

    // Calcul de la fausse diagonale : on ne passera la que si min et max sont différents.
    for(int i = min; i < max; i++){
        double curr_dist = 0;
        if(n1 < n2){
            curr_dist = dist_e(traj1+2*(min - 1), traj2+2*i);
        }else{
            curr_dist = dist_e(traj1+2*i, traj2+2*(min - 1));
        }
        distances[(i * min) + (min - 1)] = curr_dist;

        if(curr_dist > max_diag){
            max_diag = curr_dist;
        }
    }
    
    // Parcours triangulaire superieur.
    for(int i = 0; i < n1; i++){
        bool bigger = false;
        for(int j = i + 1; j < n2; j++){
            if(bigger){
                distances[i * n2 + j] = INFINITY;
            }else{
                double curr_dist = dist_e(traj1+2*i, traj2+2*j);
                if(curr_dist <= max_diag){
                    distances[i * n2 + j] = curr_dist;
                }else{
                    bigger = true;
                    distances[i * n2 + j] = INFINITY;
                }
            }
        }
    }

    // Parcours triangulaire inférieur
    for(int i = 0; i < n2; i++){
        bool bigger = false;
        for(int j = i; j < n1; j++){
            if(bigger){
                distances[j * n2 + i] = INFINITY;
            }else{
                double curr_dist = dist_e(traj1+2*j, traj2+2*i);
                if(curr_dist <= max_diag){
                    distances[j * n2 + i] = curr_dist;
                }else{
                    bigger = true;
                    distances[j * n2 + i] = INFINITY;
                }
            }
        }
    } */

    // Matrice de distance non optimisée

   /* for (int i=0; i<n1; i++){
        for (int j=0; j<n2; j++){
            double curr_dist = dist_e(traj1+2*i, traj2+2*j);
            distances[i * n2 + j] = curr_dist;
        }
    }*/

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
        traj1 = lecture(argv[1],&n_traj[0]);
        traj2 = lecture(argv[2],&n_traj[1]);

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

    MPI_Bcast(n_traj, 2, MPI_INT, 0, MPI_COMM_WORLD);

    if(pid != 0){
        traj1 = new double[2 * n_traj[0]];
    }

    MPI_Bcast(traj1, 2 * n_traj[0], MPI_DOUBLE, 0, MPI_COMM_WORLD);

    int nlocal = n_traj[1] / nprocs;
    int reste = n_traj[1] % nprocs;

    if(pid < reste){
        nlocal ++;
    }

    nlocal *= 2;

    traj2_scattered = new double[nlocal];

    MPI_Scatterv(traj2, sendcounts, displs, MPI_DOUBLE, traj2_scattered, nlocal, MPI_DOUBLE, 0, MPI_COMM_WORLD);

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

    MPI_Gatherv(mat_dist, (n_traj[0] * (nlocal / 2)) , MPI_DOUBLE, distances, distances_recvcounts, distances_displs, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    if(pid == 0){
        // Transposition du tableau 
        for(int i = 0; i < n_traj[0]; i++){
            for(int j = i+1; j < n_traj[1]; j++){
                double temp = distances[i * n_traj[0] + j];
                std::cout << i * n_traj[0] + j << " avec " << j * (n_traj[1] + 1) + i << std::endl;
                distances[i * n_traj[0] + j] = distances[i * n_traj[0] + j + i];
                distances[j * (n_traj[1] + 1) + i] = temp;
            }
        }
        for(unsigned int i = 0; i < n_traj[0] * n_traj[1]; i++){
            std::cout << distances[i] << "   ";
            if(i % n_traj[1] == 0  && i != 0){
                std::cout << std::endl;
            }
        }
    }
    /*if(pid == 0){
        cout << "la matrice de distance" << endl;

        for (int i=0; i<n_traj1; i++) {
            for (int j = 0; j < n_traj2; j++)
                cout << mat_dist[i * n_traj2 + j] << "  ";
            cout << endl;
        }
    }


    double* frechet = matrice_frechet(mat_dist,n_traj1,n_traj2);

    cout << "la matrice de fréchet" << endl;

    for (int i=0; i<n_traj1; i++) {
        for (int j = 0; j < n_traj2; j++)
            cout << frechet[i * n_traj2 + j] << " ";
        cout << endl;
    }
    */ 
    MPI_Finalize();

   //delete[] traj2_scattered;

    if(pid == 0){
        delete[] displs;
        delete[] sendcounts;
        delete[] distances;
        delete[] distances_displs;
    }

    delete[] mat_dist;

    return 0;
}
