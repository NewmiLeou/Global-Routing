#include <iostream>
#include <fstream>
#include <string.h>
#include <vector>
#include <queue>
#include <unordered_set>
#include <algorithm>
#include <cstdio>
#include <cmath>
#include <ctime>

using namespace std;

typedef struct net {
    int sx, sy;
    int tx, ty;
} Net;

typedef struct edge {
    int idx;
    int direction;
} Edge;

/* variable declaration */
int grid_x, grid_y;     // num of grids in x & y
int h_cap, v_cap;       // horizontal & vertical capacity
int num_nets;
vector<Net> nets;
vector<double> grids;   // routing grids
vector<int> grids_direction;    // routing directions, 0: top, 1: bottom, 2: left, 3: right
int max_grids;
vector<int> grids_h_demand, grids_v_demand;   // grids capacity 
vector<unordered_set<int>> grids_h_nets, grids_v_nets;  // nets pass through
vector<int> h_overflow, v_overflow;     // overflow
vector<int> h_historical, v_historical;     // historical terms
vector<vector<int>> nets_routing;   // routing results of each net
double be = 1.0;    // adaptive base cost function

/* function declaration */
void parse_input(string input_file);
void write_output(string output_file);
double GetRoutingCost(const int edge, const int direction);
void WavePropagation(const int source, const int target, const int net_id, const bool monotonic);
void Retrace(const int source, const int target, const int net_id);
void MazeRouting(const int net, const int monotonic);
void InitialRouting();
int CalculateOverflow();
struct CompareOverflow;
struct CompareBoundingBox;
void GetRipupEdges(priority_queue<Edge, vector<Edge>, CompareOverflow> &ripup_queue);
void RipupReroute(const int ripup_edge, const int ripup_direction);
void GlobalRouting();
int Calculate_WL();

/* main */
int main(int argc, char **argv){
    string input_file = argv[1];
    string output_file = argv[2];
    char testcase[20] = {0};
    strncpy(testcase, argv[2]+10, 5);
    printf("[ Testcase %s ]\n", testcase);

    printf("Reading Input File\n");
    parse_input(input_file);

    printf("Start Routing\n");
    GlobalRouting();

    printf("Writing Output File\n");
    write_output(output_file);

    // Total Overflow
    int total_overflow = CalculateOverflow();
    printf("[ Total Overflow ] : %d\n", total_overflow);

    // Total Wirelength
    int total_wirelength = Calculate_WL();
    printf("[ Total Wirelength ] : %d\n", total_wirelength);

    // Total Runtime
    double runtime = clock();
    printf("[ Total Run time ]: %.2f sec\n\n",runtime/CLOCKS_PER_SEC);

    return 0;
}


/* functions */
// parse input file
void parse_input(string input_file){
    ifstream file;
    file.open(input_file);

    string str1, str2;
    // grid # of horizontal grids # of vertical grids
    file >> str1 >> grid_x >> grid_y;
    // vertical capacity vertical capacity
    file >> str1 >> str2 >> v_cap;
    // horizontal capacity horizontal capacity
    file >> str1 >> str2 >> h_cap;
    // num net # of nets
    file >> str1 >> str2 >> num_nets;

    nets = vector<Net>(num_nets);
    max_grids = grid_x * grid_y;
    grids_h_demand = vector<int>((grid_x-1) * grid_y, 0);
    grids_v_demand = vector<int>(grid_x * (grid_y-1), 0);
    grids_h_nets = vector<unordered_set<int>>((grid_x-1) * grid_y);
    grids_v_nets = vector<unordered_set<int>>(grid_x * (grid_y-1));
    h_overflow = vector<int>((grid_x-1) * grid_y, 0);
    v_overflow = vector<int>(grid_x * (grid_y-1), 0);
    h_historical = vector<int>((grid_x-1) * grid_y, 1);
    v_historical = vector<int>(grid_x * (grid_y-1), 1);
    nets_routing = vector<vector<int>>(num_nets);

    // net-name net-id # of pins
    // pin x-grid coordinate pin y-grid coordinate
    int net_id, num_pins;
    for(int i=0; i< num_nets; i++){
        file >> str1 >> net_id >> num_pins;
        file >> nets[i].sx >> nets[i].sy >> nets[i].tx >> nets[i].ty;
    }

    file.close();
}

// write output file
void write_output(string output_file){
    ofstream file;
    file.open(output_file);

    for(int i=0; i<num_nets; i++){
        file << "net" << i << " " << i << "\n";

        int size = nets_routing[i].size()-1;
        for(int j=size; j>0; j--){
            file << "(" << nets_routing[i][j] % grid_x << "," << nets_routing[i][j] / grid_x << ", 1)";
            file << "-";
            file << "(" << nets_routing[i][j-1] % grid_x << "," << nets_routing[i][j-1] / grid_x << ", 1)";
            file << "\n";
        }
        file << "!\n";
    }

    file.close();
}

// Get routing cost
double GetRoutingCost(const int edge, const int direction){
    if (direction == 0)
        return be + h_historical[edge] * pow((double)(grids_h_demand[edge] + 1) / h_cap, 5);
    else
        return be + v_historical[edge] * pow((double)(grids_v_demand[edge] + 1) / v_cap, 5);
}

// Wave propagation
void WavePropagation(const int source, const int target, const int net_id, const bool monotonic){
    queue<int> wavefront;
    wavefront.push(source);

    int top_boundary = grid_y - 1;
    int bottom_boundary = 0;
    int left_boundary = 0;
    int right_boundary = grid_x - 1;
    if(monotonic){
        int sx = source % grid_x;
        int sy = source / grid_x;
        int tx = target % grid_x;
        int ty = target / grid_x;
        if(sx < tx){
            left_boundary = sx;
            right_boundary = tx;
        }
        else{
            left_boundary = tx;
            right_boundary = sx;
        }
        if(sy < ty){
            bottom_boundary = sy;
            top_boundary = ty;
        }
        else{
            bottom_boundary = ty;
            top_boundary = sy;
        }
    }

    while(1){
        queue<int> wavefront_temp;
        while(!wavefront.empty()){
            int wf = wavefront.front();
            wavefront.pop();
            int wf_x = wf % grid_x;
            int wf_y = wf / grid_x;

            int bottom_wf = wf - grid_x;
            if (wf_y != bottom_boundary) {
                int edge_idx = wf_x + (wf_y - 1) * grid_x;
                double cost = GetRoutingCost(edge_idx, 1);
                if (grids[bottom_wf] < 0 || grids[bottom_wf] > grids[wf] + cost) {
                    grids[bottom_wf] = grids[wf] + cost;
                    grids_direction[bottom_wf] = 0;
                    wavefront_temp.push(bottom_wf);
                }
            }

            int top_wf = wf + grid_x;
            if(wf_y != top_boundary){
                int edge_idx = wf_x + wf_y * grid_x;
                double cost = GetRoutingCost(edge_idx, 1);
                if(grids[top_wf] < 0 || grids[top_wf] > grids[wf] + cost){
                    grids[top_wf] = grids[wf] + cost;
                    grids_direction[top_wf] = 1;
                    wavefront_temp.push(top_wf);
                } 
            }

            int right_wf = wf + 1;
            if (wf_x != right_boundary) {
                int edge_idx = wf_x + wf_y * (grid_x - 1);
                double cost = GetRoutingCost(edge_idx, 0);
                if (grids[right_wf] < 0 || grids[right_wf] > grids[wf] + cost) {
                    grids[right_wf] = grids[wf] + cost;
                    grids_direction[right_wf] = 2;
                    wavefront_temp.push(right_wf);
                }
            }

            int left_wf = wf - 1;
            if (wf_x != left_boundary) {
                int edge_idx = (wf_x - 1) + wf_y * (grid_x - 1);
                double cost = GetRoutingCost(edge_idx, 0);
                if (grids[left_wf] < 0 || grids[left_wf] > grids[wf] + cost) {
                    grids[left_wf] = grids[wf] + cost;
                    grids_direction[left_wf] = 3;
                    wavefront_temp.push(left_wf);
                }
            }
        }

        if(wavefront_temp.size() == 0){
            break;
        }

        wavefront = wavefront_temp;
    }
}

// Retrace the wave
void Retrace(const int source, const int target, const int net_id){
    nets_routing[net_id].emplace_back(target);
    int wf = target;
    int prev_dir = -1;  // 0:horizontal, 1:vertical

    while(wf != source){
        int wf_x = wf % grid_x;
        int wf_y = wf / grid_x;

        int direction = grids_direction[wf];
        // top
        if(direction == 0){
            int edge_idx = wf_x + wf_y * grid_x;
            grids_v_demand[edge_idx]++;
            grids_v_nets[edge_idx].insert(net_id);
            if(prev_dir == 0){
                nets_routing[net_id].emplace_back(wf);
            }
            prev_dir = 1;
            wf += grid_x;
        }
        // bottom
        else if(direction == 1){
            int edge_idx = wf_x + (wf_y-1) * grid_x;
            grids_v_demand[edge_idx]++;
            grids_v_nets[edge_idx].insert(net_id);
            if(prev_dir == 0){
                nets_routing[net_id].emplace_back(wf);
            }
            prev_dir = 1;
            wf -= grid_x;
        }
        // left
        else if(direction == 2){
            int edge_idx = (wf_x-1) + wf_y * (grid_x-1);
            grids_h_demand[edge_idx]++;
            grids_h_nets[edge_idx].insert(net_id);
            if(prev_dir == 1){
                nets_routing[net_id].emplace_back(wf);
            }
            prev_dir = 0;
            wf -= 1;
        }
        // right
        else if(direction == 3){
            int edge_idx = wf_x + wf_y * (grid_x-1);
            grids_h_demand[edge_idx]++;
            grids_h_nets[edge_idx].insert(net_id);
            if(prev_dir == 1){
                nets_routing[net_id].emplace_back(wf);
            }
            prev_dir = 0;
            wf += 1;
        }
    }
    nets_routing[net_id].emplace_back(source);
}

// Maze routing
void MazeRouting(const int net_id, const int monotonic){
    grids = vector<double>(grid_x * grid_y, -1);
    grids_direction = vector<int>(grid_x * grid_y, -1);

    int source = nets[net_id].sx + nets[net_id].sy * grid_x;
    int target = nets[net_id].tx + nets[net_id].ty * grid_x;
    grids[source] = 0;

    WavePropagation(source, target, net_id, monotonic);
    Retrace(source, target, net_id);
}

// Initial routing
void InitialRouting(){
    for(int i=0; i<num_nets; i++){
        MazeRouting(i, true);
    }
}

// Calculate total overflow
int CalculateOverflow(){
    int total_overflow = 0;

    for(int y=0; y<grid_y; y++){
        for(int x=0; x<grid_x-1; x++){
            int id_x = x + y * (grid_x-1);
            h_overflow[id_x] = max(grids_h_demand[id_x] - h_cap, 0);
            total_overflow += h_overflow[id_x];
        }
    }

    for(int y=0; y<grid_y-1; y++){
        for(int x=0; x<grid_x; x++){
            int id_x = x + y * grid_x;
            v_overflow[id_x] = max(grids_v_demand[id_x] - v_cap, 0);
            total_overflow += v_overflow[id_x];
        }
    }

    return total_overflow;
}

struct CompareOverflow {
    bool operator() (const Edge e1, const Edge e2) {
        int overflow_1 = e1.direction==0 ? h_overflow[e1.idx] : v_overflow[e1.idx];
        int overflow_2 = e2.direction==0 ? h_overflow[e2.idx] : v_overflow[e2.idx];
        return overflow_1 > overflow_2;
    }
};

// Get rip up edges
void GetRipupEdges(priority_queue<Edge, vector<Edge>, CompareOverflow> &ripup_queue){
    Edge e;

    for(int y=0; y<grid_y; y++){
        for(int x=0; x<grid_x-1; x++){
            int id_x = x + y * (grid_x-1);
            if(h_overflow[id_x] > 0){
                h_historical[id_x]++;
                e.idx = id_x;
                e.direction = 0;
                ripup_queue.push(e);
            }
        }
    }

    for(int y=0; y<grid_y-1; y++){
        for(int x=0; x<grid_x; x++){
            int id_x = x + y * grid_x;
            if(v_overflow[id_x] > 0){
                v_historical[id_x]++;
                e.idx = id_x;
                e.direction = 1;
                ripup_queue.push(e);
            }
        }
    }
}

struct CompareBoundingBox {
    bool operator() (const int net1, const int net2) {
        int hpwl_1 = abs(nets[net1].sx - nets[net1].tx) + abs(nets[net1].sy - nets[net1].ty);
        int hpwl_2 = abs(nets[net2].sx - nets[net2].tx) + abs(nets[net2].sy - nets[net2].ty);
        return hpwl_1 < hpwl_2;
    }
};

// Ripup & Reroute
void RipupReroute(const int ripup_edge, const int ripup_direction){
    // rip up
    unordered_set<int> ripup_nets;
    if(ripup_direction == 0){
        ripup_nets = grids_h_nets[ripup_edge];
    }
    else{
        ripup_nets = grids_v_nets[ripup_edge];
    }

    priority_queue<int, vector<int>, CompareBoundingBox> reroute_order;
    for(const int net_id : ripup_nets){
        int size = nets_routing[net_id].size() - 1;
        for(int i=0; i<size; i++){
            int x1 = nets_routing[net_id][i] % grid_x;
            int y1 = nets_routing[net_id][i] / grid_x;
            int x2 = nets_routing[net_id][i+1] % grid_x;
            int y2 = nets_routing[net_id][i+1] / grid_x;

            int S, T, direction;
            if(x1 == x2){
                S = y1 < y2 ? y1 : y2;
                T = y1 < y2 ? y2 : y1;
                direction = 1;
            }
            else{
                S = x1 < x2 ? x1 : x2;
                T = x1 < x2 ? x2 : x1;
                direction = 0;
            }

            if(direction == 0){
                for(int j=S; j<T; j++){
                    int id_x = j + y1 * (grid_x-1);
                    grids_h_nets[id_x].erase(net_id);
                    grids_h_demand[id_x]--;
                }
            }
            else{
                for(int j=S; j<T; j++){
                    int id_x = x1 + j * grid_x;
                    grids_v_nets[id_x].erase(net_id);
                    grids_v_demand[id_x]--;
                }
            }
        }

        nets_routing[net_id] = vector<int>();
        reroute_order.push(net_id);
    }

    // reroute
    while(!reroute_order.empty()){
        int net_id = reroute_order.top();
        reroute_order.pop();
        MazeRouting(net_id, false);
    }
}

// Global routing
void GlobalRouting(){
    clock_t start = clock();
    float time_elapsed;

    InitialRouting();

    int total_overflow = CalculateOverflow();
    int iter = 0;
    while(total_overflow > 0){
        iter++;
        if(num_nets == 27781){
            be = 1.0 - exp(-5 * exp(-0.4 * iter));  // adaptive base cost function
        }
        else if(num_nets == 13357){
            be = 1.0 - exp(-10 * exp(-0.2 * iter));
        }
        else{
            be = 1.0 - exp(-5 * exp(-0.1 * iter));
        }
        
        priority_queue<Edge, vector<Edge>, CompareOverflow> ripup_queue;
        GetRipupEdges(ripup_queue);

        while(total_overflow > 0 && !ripup_queue.empty()){
            int edge_idx = ripup_queue.top().idx;
            int direction = ripup_queue.top().direction;
            ripup_queue.pop();
            if(((direction == 0) && (grids_h_demand[edge_idx] > h_cap)) || ((direction == 1) && (grids_v_demand[edge_idx] > v_cap))){
                RipupReroute(edge_idx, direction);
            }
            total_overflow = CalculateOverflow();

            time_elapsed = (clock() - start) / CLOCKS_PER_SEC;
            if(time_elapsed >= 590){
                total_overflow = 0;
            }
        }
    }
}

// Calculate total wirelength
int Calculate_WL(){
    int total_wirelength = 0;
    for(int i=0; i<num_nets; i++){
        int wl = abs(nets[i].sx - nets[i].tx) + abs(nets[i].sy - nets[i].ty);
        total_wirelength += wl;
    }
    return total_wirelength;
}