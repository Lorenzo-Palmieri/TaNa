//
//  Fast+.cpp
//
//
//  Created by Lorenzo Palmieri on 05/01/17.
//
//

#include "Tana+.hpp" // Include all libraries.
using namespace std; // avoid using std:: along the code.
// random_device rd; // obtain a random number from hardware.
// mt19937 eng(rd()); // seed the generator.
mt19937 eng(12345);
// mt19937 is Mersenne Twister pseudo-random generator of 32-bit numbers.
// N.B. in this code the genome is made of a string of 0 and 1, rather than -1 and 1.

// -----------------------------------
//       INITIALIZE PARAMETERS
// -----------------------------------
const int generations = 100000, L = 15, No = 1000, K = 1, genomes = 32768;
const double theta = 0.25, pmut = 0.0125, pkill = 0.2, c = 14.3 , mu = 0.005;
const double pi = 3.1415927;
// -----------------------------------
//        DECLARE FUNCTIONS
// -----------------------------------
double gaussian(void);
double uniform(void);
int choose_in_range(int a, int b);
void initial_pop( array<int, genomes> &vec);
void Kill (array<int,genomes> &spec, vector<int> &ex, int &N);
void get_row(int (&m)[K][L], int row, array<int, L> &r);
void Baby(array<int,genomes> &spec, vector<int> &ex, double (&J)[genomes][genomes], int &N);
void compute_interactions(int x, int y, array<double,genomes> &A, array <int, genomes> &B, array<double,2> &inter);
void Bin_recursive(int n, int j, array<int, L> &vec);
void ConvertToBinary(int n, int j, array<int, L> &vec);
int ConvertToDec(array<int, L> &arr);
void Create_J(double (&J)[genomes][genomes]);
// -----------------------------------
//            TEMPLATES
// -----------------------------------
template <class T>
T zerovec (T &vec) {
    int i;
    for (i=0; i < vec.size(); i++) {
        vec[i]=0;
    }
    return (vec);
}

    static array <int, genomes> species;
    static double J[genomes][genomes]={};

// -----------------------------------
//                MAIN
// -----------------------------------
int main(){

    vector<int> existent;
    int i,N, count_gen, step;
    double tau;
    
    Create_J(J);
    // set all elements of species to zero and initialize population at random.
    zerovec(species);
    initial_pop(species);
    N = No;
    tau = round(double(N)/pkill);
    // Initialize a vector which contains the label of all existent species.
    for (i=0; i<genomes; i++) {
        if (species[i] != 0) {
            existent.push_back(i);
        }
    }
    // create output files.
    ofstream mydata;
    mydata.open("data.txt");
    ofstream spec;
    spec.open("species.txt");
    ofstream allspec;
    allspec.open("complete.txt");
    
    count_gen = 0;
    step = 0;
    while (count_gen < generations) {
        step++;
        // Try to kill an individual.
        Kill (species, existent, N);
        // Try to reproduce an individual
        Baby(species, existent, J, N);
        if(step == tau){
            count_gen++;
            step=0;
            tau = round(double(N)/pkill);
            mydata << count_gen << "\t" << N << endl;
            for (i=0; i<existent.size(); i++) {
                spec << existent[i] << "\t";
                allspec << species[existent[i]] << "\t";
            }
            spec << endl;
            allspec << endl;
            
            if (count_gen % (generations/10) == 0) {
                cout << count_gen << endl;
            }
        }
    }
    mydata.close();
    spec.close();
    allspec.close();
}

// -----------------------------------
//             FUNCTIONS
// -----------------------------------

void Baby(array<int,genomes> &spec, vector<int> &ex, double (&J)[genomes][genomes], int &N) {
    int chosen, selected, add, i, s1, s2;
    double t1=0, H, poff;
    array <int,L > b1, b2;
    
    // select an element of existent.
    selected = choose_in_range(0, ex.size() - 1);
    // associate the correspondent species.
    chosen = ex[selected];
    // Compute poff
    for (i=0; i < ex.size(); i++) {
        t1 += J[chosen][ex[i]] * spec[ex[i]];
    }
    H = c*t1/N - mu*N;
    poff = 1/(1+exp(-H));
    // try to reproduce a species.
    if (uniform() <= poff){
        //total population increases by 1.
        N++;
        // Kill the parent.
        spec[chosen] -= 1;
        // check if it implies extinction.
        if (spec[chosen] == 0) {
            ex.erase (ex.begin() + selected);
        }
        ConvertToBinary(chosen, L, b1);
        b2 = b1;
        // Apply mutations to babies.
        for (i=0; i<L; i++) {
            // 1st baby.
            if (uniform() <= pmut) {
                if (b1[i] == 1) {
                    b1[i]=0;
                } else {
                    b1[i]=1;
                }
            }
            // 2nd baby.
            if (uniform() <= pmut) {
                if (b2[i] == 1) {
                    b2[i]=0;
                } else {
                    b2[i]=1;
                }
            }
        }
        // Convert to decimal.
        s1 = ConvertToDec(b1);
        s2 = ConvertToDec(b2);
        // Check if there are new species.
        add=0;
        for (i=0; i<genomes; i++) {
            if(spec[i]==0 && (s1==i || s2==i)){
                ex.push_back(i);
                add++;
            }
            if(add==1 && s1==s2){break;}
            if(add==2){break;}
        }
        // Add babies to species.
        spec[s1] += 1;
        spec[s2] += 1;
    }
}

void Create_J(double (&J)[genomes][genomes]){
    int i,k;
    array <int, genomes> B;
    array <double, genomes> A;
    array<double,2> interactions;
    
    // Initialize vectors A and B for correlations.
    for (i=0; i < genomes; i++) {
        A[i]= gaussian();
        if(uniform() <= theta){
            B[i] = 1;
        }else{
            B[i] = 0;
        }
    }
    // create the J matrix
    for (i=1; i<genomes; i++) {
        for (k=0; k<i; k++) { // skip i=j (already zero).
            compute_interactions(i,k,A,B,interactions);
            J[i][k] = interactions[0];
            J[k][i] = interactions[1];
        }
    }
}


void compute_interactions(int x, int y, array<double,genomes> &A, array <int, genomes> &B, array<double,2> &inter ){
    double w1=0, w2=0, w3=0;
    int z,i,j;
    int z1[K][L]={},  z2[K][L]= {}, z3[K][L]= {};
    int n1[K], n2[K], n3[K];
    array<int, L> v1,v2,v3,r;
    
    z=x^y; // XOR
    // convert to binary.
    ConvertToBinary(x, L, v1);
    ConvertToBinary(y, L, v2);
    ConvertToBinary(z, L, v3);
    // split to create correlations.
    for (j=0; j<K; j++) {
        for (i=0; i<L/K; i++) {
            z1[j][i + j*L/K] = v1[i + j*L/K];
            z2[j][i + j*L/K] = v2[i + j*L/K];
            z3[j][i + j*L/K] = v3[i + j*L/K];
        }
        // convert to decimal.
        get_row(z1,j,r);
        n1[j] = ConvertToDec(r);
        get_row(z2,j,r);
        n2[j] = ConvertToDec(r);
        get_row(z3,j,r);
        n3[j] = ConvertToDec(r);
        
    }
    if (B[z]==0 || x==y) {
        inter[0]=0;
        inter[1]=0;
    } else{
        for (j=0; j<K; j++) {
            w1 +=  A[n1[j]]/sqrt(double(K));
            w2 +=  A[n2[j]]/sqrt(double(K));
            w3 +=  A[n3[j]]/sqrt(double(K));
        }
        inter[0] = w1 * w3;
        inter[1] = w2 * w3;
    }
}

void get_row(int (&m)[K][L], int row, array<int, L> &r) {
    int i;
    for (i=0; i<L; i++) {
        r[i] = m[row][i];
    }
}

void Kill (array<int,genomes> &spec, vector<int> &ex, int &N){
    int chosen, selected;
    // select an element of existent.
    selected = choose_in_range(0, ex.size() - 1);
    // associate the correspondent species.
    chosen = ex[selected];
    // Attempt to kill an individual.
    if (uniform() <= pkill) {
        // total population decreases by 1.
        N--;
        spec[chosen] -= 1;
        // Is the chosen species extinct?
        if (spec[chosen] == 0) {
            ex.erase(ex.begin() + selected);  // being "ex" small, is better than using std::list.
        }
    }
    // check if all species are dead.
    if (ex.size() == 0) {
        cout << "extinction" << endl;
        exit(EXIT_FAILURE);
    }
}


void initial_pop( array<int, genomes> &vec){
    int i, select;
    for (i=1; i <= No; i++) {
        select = choose_in_range(0,genomes-1);
        vec[select] += 1;
    }
}

double gaussian(void){
    double u1=0, u2=0;
    double gauss_var;
    // exclude the zero values.
    while (u1 == 0 || u2 == 0) {
        u1 = uniform();
        u2 = uniform();
    }
    gauss_var = sqrt( -2.0*log(u1)) * cos(2.0*pi*u2);
    return(gauss_var);
}


double uniform(void){
    uniform_real_distribution<> rand(0, 1); // define the range [a,b], extremes included.
    return rand(eng);
}


int choose_in_range(int a, int b){
    uniform_int_distribution<> choose(a, b); // define the range [a,b], extremes included.
    return choose(eng);
}


void Bin_recursive(int n, int j, array<int, L> &vec){ // It must be called setting j=L.
    j--;
    if (n / 2 != 0) {
        Bin_recursive(n/2, j, vec);
    }
    vec[j] = n % 2;
}


void ConvertToBinary(int n, int j, array<int, L> &vec){
    int i;
    if(n >= pow(2,L)){
        cout << "the number cannot be represented" << endl;
        exit (EXIT_FAILURE);
    }
    for (i=0; i<L; i++) {
        vec[i]=0;
    }
    Bin_recursive(n, j, vec);
}


int ConvertToDec(array<int, L> &arr){
    int i, dec=0;
    for (i=0; i<L;i++) {
        if(arr[L-1-i] == 1){
            dec += pow(2,i);
        }
    }
    return dec;
}
