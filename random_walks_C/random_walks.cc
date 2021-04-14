#include <cmath>
#include <iostream>
#include <vector>
#include <time.h>
#include <chrono>
#include <fstream>
#include <cstdlib>
#include <string>

using namespace std;

// Program parameters: iterations, n_walkers, rec_lvl

// (1, 10^4), (2, 10^5), (3, 10^6), (4, 10^7), (5,10^7), (6, 10^8), (7, 2*10^8)

int iterations = 2 * pow(10, 8);
int n_walkers = 100;
float data_size = 10000;
int rec_lvl = 7;
float L = 1;
float Rmax = L / 2;
float N = pow(3, rec_lvl);
float tam_min = L / pow(3, rec_lvl); 
float dt = pow(tam_min, 2) / 10;
float sqrt_dt = sqrt(dt);

float default_filename = false;
string data_folder = "";
vector<vector<int>> grid;

struct Walker {
    float x, y;
    float x_pbc, y_pbc;
    float x_init, y_init;
};

int isSierpinskiCarpetPixelFilled(int x, int y, int width, int height, int d, int max_depth, float v) {
	int x2 = x * 3 / width;
	int y2 = y * 3 / height;
	if (x2 == 1 && y2 == 1) {
        return v;
    } 

	// if ((x <= 0)||(y <= 0)||(x>=width)||(y>=height)) 
	// 	return 0;

    if (d < max_depth) {
        x -= (x2 * width + 2) / 3; 
        y -= (y2 * height + 2) / 3;
        width = (width + 2 - x2) / 3;
        height = (height + 2 - y2) / 3;
	    return isSierpinskiCarpetPixelFilled(x, y, width, height, d + 1, max_depth, v);
    }
    else {
        return 0;
    }
}


void sierpinski(vector<vector<int>> &mat, int rec, bool verbose){
    int rows = mat.size();
    int cols = mat[0].size();
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            mat[i][j] = isSierpinskiCarpetPixelFilled(i, j, rows, cols, 1, rec, 1);
        }
        if (verbose && i % 500 == 0) cout << "fila " << i << endl;
    }
}


void print_mat(vector<vector<int>> &mat){
    int rows = mat.size();
    int cols = mat[0].size();
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            if (mat[i][j]) cout << 'x';
            else cout << '.';
            cout << ' ';
        }
        cout << endl;
    }
}


bool valid_pos(float x, float y) {
    int g_posx = (x + Rmax) * N/L;
    int g_posy = (y + Rmax) * N/L;
    if (g_posy < 0 || g_posx < 0 || g_posy >= N || g_posx >= N)
        cout << "EXCEPTION(Valid pos): " << x << ' ' << y << ' ' << g_posx << ' ' << g_posy << endl;
    bool res = grid[g_posx][g_posy] == 0; 
    return res;
}


void init_walkers(vector<Walker> &walkers) {
    for (int i = 0; i < walkers.size(); ++i) {
        float x = float(rand())/(RAND_MAX + 1) * L - Rmax;
        float y = float(rand())/(RAND_MAX + 1) * L - Rmax;
        
        while (!valid_pos(x, y)) {
            x = float(rand())/RAND_MAX * L - Rmax;
            y = float(rand())/RAND_MAX * L - Rmax;
        }
        
        walkers[i].x = x;
        walkers[i].x_pbc = x;
        walkers[i].y = y;
        walkers[i].y_pbc = y;
        walkers[i].x_init = x;
        walkers[i].y_init = y;
    }
}


void random_values(float &z1, float &z2) {
    float phi = 2 * M_PI * float(rand())/RAND_MAX;
    float rand_val = float(rand())/RAND_MAX;
    if (rand_val == 0) rand_val = 1.0/RAND_MAX;
    float r = sqrt(-2 * log(rand_val));
    z1 = r * cos(phi);
    z2 = r * sin(phi);
    if (r > 10000 || r < -10000 || phi > 10000 || phi < -10000)
        cout << "rand: "<< r  << ' ' << rand_val << ' ' << phi << endl;
}


void update_walker(Walker *walker) {
    float a, b;
    random_values(a, b);

    float tmp_x = walker->x + a * sqrt_dt;
    float tmp_x_pbc = walker->x_pbc + a * sqrt_dt;
    float tmp_y = walker->y + b * sqrt_dt;
    float tmp_y_pbc = walker->y_pbc + b * sqrt_dt;

    int g_posx = (tmp_x_pbc + Rmax) * N/L;
    int g_posy = (tmp_y_pbc + Rmax) * N/L;

    if (g_posx < 0) {
        tmp_x_pbc += L; 
        g_posx += N;
    }
    else if (g_posx >= N) {
        tmp_x_pbc -= L; 
        g_posx -= N;
    }
    
    if (g_posy < 0) {
        tmp_y_pbc += L; 
        g_posy += N;
    } 
    else if (g_posy >= N) {
        tmp_y_pbc -= L; 
        g_posy -= N;
    }

    // valid_pos(tmp_x_pbc, tmp_y_pbc) 
    if (grid[g_posx][g_posy] == 0) {
        walker->x = tmp_x;
        walker->y = tmp_y;
        walker->x_pbc = tmp_x_pbc;
        walker->y_pbc = tmp_y_pbc;
    }
}


void storing_distances(const vector <float> &distances) {
    string file_name = data_folder + '/';
    if (default_filename) 
        file_name += "distances.txt";
    else
        file_name += "distances_" + to_string(int(rec_lvl)) + '_' + to_string(int(n_walkers)) + '_' + to_string(int(N)) + ".txt";

    cout << "Storing distances on " << file_name << endl;
    ofstream dist_file(file_name);

    dist_file << iterations << ' ' << n_walkers << ' ' << rec_lvl << ' ';
    dist_file << L << ' ' << N << ' ' << tam_min << ' ' << dt << ' ' << data_size << endl;  

    dist_file << distances[0] / n_walkers;
    int batch_size = iterations / data_size;

    for (int i = 0; i < distances.size(); i += batch_size) {
        dist_file << ',' << distances[i] / n_walkers;
    }

    dist_file.close();
}


void storing_track_mat(const vector<vector<int>> &track_mat) {
    string file_name = data_folder + '/';
    if (default_filename)
        file_name += "track_mat.txt";
    else
        file_name += "track_mat_" + to_string(int(rec_lvl)) + '_' + to_string(int(n_walkers)) + '_' + to_string(int(N)) + ".txt";
        
    cout << "Storing track_mat on " << file_name << endl;

    ofstream track_mat_file(file_name);

    int rows = track_mat.size();
    int cols = track_mat[0].size();
    for (int i = 0; i < rows; ++i) {
        track_mat_file << track_mat[i][0];
        for (int j = 1; j < cols; ++j) {
            track_mat_file << ',' << track_mat[i][j];
        }
        track_mat_file << endl;
    }

    track_mat_file.close();
}


int main(int argc, char* argv[]) {
    // Read parameters
    if (argc > 1) {
        iterations = atoi(argv[1]);
        n_walkers = atoi(argv[2]);
        rec_lvl = atoi(argv[3]);
        data_folder = argv[4];

        // Compute new global variables        
        N = pow(3, rec_lvl);
        tam_min = L / pow(3, rec_lvl); 
        dt = pow(tam_min, 2) / 10;
        sqrt_dt = sqrt(dt);
        cout << "Parameters: " << iterations << ' ' << n_walkers << ' ' << rec_lvl << endl;
    }

    data_size = data_size > iterations ? iterations : data_size;


    using std::chrono::high_resolution_clock;
    using std::chrono::duration_cast;
    using std::chrono::duration;
    using std::chrono::milliseconds;
    auto t1 = high_resolution_clock::now();

    srand(time(0));

    grid = vector<vector<int>> (N, vector<int> (N, 0));
    vector<vector<int>> track_mat(N, vector<int> (N, 0));
    vector <float> distances (iterations, 0);
    
    cout << "Computing sierpinski..." << endl;
    sierpinski(grid, rec_lvl, false);

    cout << "Initialize walkers..." << endl;
    vector <Walker> walkers (n_walkers);
    init_walkers(walkers);


    auto t2 = high_resolution_clock::now();
    duration<double, std::milli> ms_double = t2 - t1;
    std::cout << "Initialization time: " << ms_double.count()/1000 << "s" << endl;

    for (int i = 0; i < n_walkers; ++i) {
        // if (i % 10 == 0) {
        //     auto elapsed = high_resolution_clock::now();
        //     ms_double = elapsed - t2;
        //     cout << "Walker: " << i << " Elapsed time: " << ms_double.count() << "ms" << endl; 
        // } 
        
        for (int j = 0; j < iterations; ++j) {
            // if (j % 50000 == 0) cout << "Iteration: " << j << endl; 
            update_walker(&walkers[i]);
            
            float dx = walkers[i].x_init - walkers[i].x;
            float dy = walkers[i].y_init - walkers[i].y;

            int g_posx = (walkers[i].x_pbc + Rmax) * N/L;
            int g_posy = (walkers[i].y_pbc + Rmax) * N/L;

            if (g_posy < 0 || g_posx < 0 || g_posy >= N || g_posx >= N)
                cout << "Update: " << walkers[i].x_pbc << ' ' << walkers[i].y_pbc << ' ' << g_posx << ' ' << g_posy << endl;
            
            track_mat[g_posx][g_posy] += 1;
            distances[j] += dx*dx + dy*dy;

        }
    }
    
    auto t3 = high_resolution_clock::now();
    ms_double = t3 - t2;
    std::cout << "Iterations time: " << ms_double.count()/1000 << "s" << endl;

    storing_distances(distances);

    storing_track_mat(track_mat);

    auto t4 = high_resolution_clock::now();
    ms_double = t4 - t3;
    std::cout << "Storing time: " << ms_double.count()/1000 << "s" << endl;

}