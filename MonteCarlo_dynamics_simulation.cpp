#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include <iterator>
#include <string>
#include <string.h>
#include <random>
#include<cmath>
#include <complex>
#include <stdlib.h>
#include <thread>
#include "func.h"


using namespace std;

enum room_size{ KIDS=2500, M=1, K_B=1};

const int STEPS = 5;
const int UNIT_TIME=1000;
const int MAX_SWAPS=10000000;

const double D=0.2;
const double EPSI=0.2;

const double A=1.1;
const double B = 1.1*0.7; // initial value of A*R0
const double R0=B/A;
const double Rc = 1/A;
const double SIGMA=0.5;
const double N=2;
const int CUT = 2000000;

const double SIZE_FACTOR = M_PI*SIGMA*SIGMA*2500;
static double LENGTH;
static double WIDTH;
static double DENSITY;


typedef complex<double> complexd;

class Variable
{
private:
	double val;
	const char* type;
public:
	double get_val() { return val; };
	Variable()
	{
		std::cout << "Create an empty variale!" << std::endl;
	};
	Variable(double val1, const char* type1)
	{
		val = periodic_boundary(val1, type1);
	}
	double periodic_boundary(double var_val, const char* unknow_type)
	{
		if (!strcmp(unknow_type, "x")) check1D(var_val, LENGTH, 0);
		if (!strcmp(unknow_type, "y")) check1D(var_val, WIDTH, 0);
		if (!strcmp(unknow_type, "a")) check1D(var_val, M_PI, -M_PI);
		return var_val;
	}

	void check1D(double &var, double upper_bound, double lower_bound)
	{
		var -= lower_bound;
		var = fmod(var, upper_bound - lower_bound);
		var += lower_bound;
		if (var < lower_bound) var += (upper_bound - lower_bound);
	}
	Variable operator+(const Variable &variable2)
	{
		Variable variable(this->val + variable2.val, this->type);
		return variable;
	};
};

struct User
{
	double x, y, a;
};

vector<User> InitUsers()
{
	ifstream file("./ini_pos.dat");
	vector<User> users;
	string line;
	User user;
	for (int i=0; i<KIDS; ++i)
	{
		getline(file, line);
		istringstream is(line);
		if (line.size())
		{
			vector<std::string> tokens;
			copy(istream_iterator<string>(is),
		    istream_iterator<string>(),
		    back_inserter(tokens));
		    Variable var_x(stod(tokens[0]), "x"), var_y(stod(tokens[1]), "y"), var_a(stod(tokens[2]), "a");
		    user.x = var_x.get_val();
		    user.y = var_y.get_val();
		    user.a = var_a.get_val();

		    users.push_back(user);
		}
	}
	return users;
};

double ModelPotential(vector2& distances)
{
	return - D * (exp(-2 * A * (distances.val1 - R0)) - 2 * exp(-A * (distances.val1 - R0))) * cos(distances.val2);
};

double PotentialForce(vector2& distances)
{
	return 2.0 * A * D * ( exp(-A*(distances.val1-R0)) - exp(-2*A*(distances.val1-R0)) ) * cos(distances.val2);
};

double SoftForce(double ds)
{
	if (ds <= SIGMA) { return N*EPSI*pow(SIGMA/ds, N); } else { return 0; }
};

double SoftPotential(double ds)
{
	if (ds <= SIGMA){ return EPSI * pow(SIGMA / ds, N); } else { return 0; };
};

void clct_dis(User& user1, User& user2, vector2& distances)
{
	distances.val1 = sqrt((user1.x-user2.x)*(user1.x-user2.x) + (user1.y-user2.y)*(user1.y-user2.y));
	distances.val2 = user1.a - user2.a;
};

double TotalPotential(vector<User>& users)
{
	double potential = 0;
	double pair_potential; vector2 distances;
	for (int i = 0; i < KIDS; i ++)
	{
		for (int j = i+1; j < KIDS; j ++)
		{
			User user_i = users[i];
			User user_j = users[j];
			clct_dis(user_i, user_j, distances);
			pair_potential = ModelPotential(distances) + SoftPotential(distances.val1);
			potential += pair_potential;
		}
	}
	return potential;
};

double TotalPressure(vector<User>& users)
{
	double P = 0;
	double pair_P; vector2 distances;
	for (int i = 0; i < KIDS; i ++)
	{
		for (int j = i+1; j < KIDS; j ++)
		{
			User user_i = users[i];
			User user_j = users[j];
			clct_dis(user_i, user_j, distances);
			pair_P = PotentialForce(distances) + SoftForce(distances.val1);
			P += pair_P;
		}
	}
	return P/(2*LENGTH*WIDTH);
};


vector2 ClctStatus(vector<User>& users, int rand_kid_idx, User& random_mover)
{
	double potential = 0;
	double pressure = 0;
	double pair_potential, pair_pressure; vector2 distances;
	for (int i = 0; i < KIDS; ++i)
	{
		if (i != rand_kid_idx)
		{
			User other_kid = users[i];
			clct_dis(random_mover, other_kid, distances);
			pair_potential = ModelPotential(distances) + SoftPotential(distances.val1);
			pair_pressure = (PotentialForce(distances) + SoftForce(distances.val1))*distances.val1;
			potential += pair_potential;
			pressure += pair_pressure;
		}
	}
	vector2 status;
	status.val1 = potential;
	status.val2 = pressure;
	return status;
};


bool UpdateStatus(double dV, double T)
{
	if (dV < 0) { return true; }
	float p;
	p = r4_uniform_ab(0, 1, seed);
	if (exp(-dV/T) > p) {
		return true;
	} else { return false; }
};

struct Energy
{
	double E, totE, sum_varE;
};

struct Count
{
	int accept_count, tot, tot_steps;
};

struct Jump
{
	double dX, dY, dA;
};

class t_star
{
	complexd fourierCos;
	complexd fourierSin;
	double tot;
public:
	void init()
	{
		fourierCos.real(0);
		fourierCos.imag(0);
		fourierSin.real(0);
		fourierSin.imag(0);
		tot = 0;
	}
	void add(vector<User>& users, double kx, double ky)
	{
		double x, y, theta, kr;
		for (int i = 0; i < users.size(); ++i)
		{
			x = users[i].x;
			y = users[i].y;
			theta = users[i].a;
			kr = kx*x + ky*y;
			complexd expTerm(cos(kr), -sin(kr));
			fourierCos += cos(theta)*expTerm;
			fourierSin += sin(theta)*expTerm;
			tot ++;
		}
	}
	vector2 findM(vector2& fourierVals)
	{
		fourierCos /= tot;
		fourierSin /= tot;

		fourierVals.val1 = fourierCos.real()*fourierCos.real() + fourierCos.imag()*fourierCos.imag();
		fourierVals.val2 = fourierSin.real()*fourierSin.real() + fourierSin.imag()*fourierSin.imag();
		return fourierVals;
	}
};

double findMmax(vector<User>& users)
{
	int max_n = 10; double kx, ky; double max_m2 = 0;
	t_star star; vector2 fourierVals;
	for (int i = -10; i <= 2*max_n; ++i)
	{
		if (i != 0) kx = 2*M_PI*i/LENGTH; else kx = 0;
		for (int j = 0; j < 2*max_n; ++j)
		{
			if (j != 0) ky = 2*M_PI*j/WIDTH; else ky = 0;
			star.init();
			star.add(users, kx, ky);
			vector2 fourierVals = star.findM(fourierVals);
			if (max_m2 < fourierVals.val1 + fourierVals.val2)
			{
				max_m2 = fourierVals.val1 + fourierVals.val2;
			}
		}
	}
	return max_m2;
}

void output(const vector<User>& users, double T)
{
	string T_str;
	stringstream num;
	num << T;
	num >> T_str;
	ofstream save1("status_" + T_str +".dat");
	for (int kid = 0; kid < KIDS; ++kid)
	{
		save1 << users[kid].x << "\t" << users[kid].y << "\t" << users[kid].a << endl;
	}
};

void update(vector<User>& users, double T, Energy& energy, vector2& P, Count& cnt, Jump &jmp, int k)
{
	double dx, dy, da, dp, p0, dE, p, heatcapacity, accept_rate, max_m2;
	int rand_kid_idx; vector2 curr_status, possible_status;
	User random_mover;
	p0 = K_B*T*KIDS/(LENGTH*WIDTH);
    string T_str; stringstream num;
	num << T; num >> T_str;
	ofstream save("output_" + T_str + ".csv", ios_base::app);

	for (int i=0; i<STEPS; i++)
	{
		cnt.tot_steps ++;

		rand_kid_idx = i4_uniform_ab(0, KIDS-1, seed); // sample a random kid using random distribution
		curr_status = ClctStatus(users, rand_kid_idx, users[rand_kid_idx]);

		dx = r4_uniform_ab(-jmp.dX, jmp.dX, seed), dy = r4_uniform_ab(-jmp.dY, jmp.dY, seed), da = r4_uniform_ab(-jmp.dA*M_PI/180, jmp.dA*M_PI/180, seed);
		Variable update_x(users[rand_kid_idx].x+dx, "x"), update_y(users[rand_kid_idx].y+dy, "y"), update_a(users[rand_kid_idx].a+da, "a");
		random_mover.x = update_x.get_val();
		random_mover.y = update_y.get_val();
		random_mover.a = update_a.get_val();

		possible_status = ClctStatus(users, rand_kid_idx, random_mover);

		dE = possible_status.val1 - curr_status.val1;
		if (UpdateStatus(dE, T))
		{
			cnt.accept_count ++;
			energy.E += dE;
			dp = (possible_status.val2 - curr_status.val2)/(2*LENGTH*WIDTH);
			P.val2 += dp;

			users[rand_kid_idx] = random_mover;
		}
	if (cnt.tot_steps%20000 == 0 && 1.0 * cnt.accept_count / cnt.tot_steps < 0.3) jmp.dX *= 0.9; jmp.dY *= 0.9; jmp.dA *= 0.9;
        if (k < CUT) continue;

        if (cnt.tot_steps%1000000 == 0) output(users, T);

        cnt.tot ++;
	    energy.totE += energy.E; energy.sum_varE += (energy.E-energy.totE/cnt.tot)*(energy.E-energy.totE/cnt.tot);
	    P.val1 += P.val2;

        if (cnt.tot_steps%UNIT_TIME == 0)
		{
			max_m2 = findMmax(users);
			p = p0 + P.val1/cnt.tot;
			heatcapacity = energy.sum_varE/(K_B*K_B*T*T)/cnt.tot;
			if (heatcapacity < 0)
			{
				cout << "err--> negative heatcapacity " << i << '\t' << heatcapacity << '\t' << energy.totE << '\t' << energy.sum_varE << '\t' << cnt.tot << endl;
			}
			accept_rate = 1.0 * cnt.accept_count / cnt.tot_steps;
            if (accept_rate < 0.3) jmp.dX *= 0.9; jmp.dY *= 0.9; jmp.dA *= 0.9;

			save << energy.E << ',' << heatcapacity/KIDS << ',' << max_m2 << ','<< p << ',' << accept_rate << endl;
		}
	}
}

bool SwapProb(double dbeta, double dE)
{
	double p = 1; double swap_p = exp(dbeta*dE);
	if (swap_p < 1)  p = swap_p;
	float rand_p = r4_uniform_ab(0, 1, seed);
	if (rand_p < p) return true;
	return false;
}

vector<double> InitTemperatures()
{
	vector<double> Temperatures;
	ifstream infile("./Temperatures.dat");
	string line;
	while (getline(infile, line))
	{
		double T = stod(line);
		Temperatures.push_back(T);
	}
	return Temperatures;
};

int main(int argc, char* argv[])
{
	DENSITY=atof(argv[1]);
	LENGTH = sqrt(SIZE_FACTOR/DENSITY);
    WIDTH = LENGTH;

	vector<User> users = InitUsers();
	double E = TotalPotential(users); double totE = 0; double sum_varE = 0;
	vector2 P; P.val1 = 0; P.val2 = TotalPressure(users); double pv; // val1 represents total pv, val2 means pv

	vector<double> Temperatures = InitTemperatures();
	cout << Temperatures.size() << endl;
	int num_ensembles = Temperatures.size();

	vector<thread> threads;
	vector<vector<User> > states(num_ensembles);
	vector<Energy> energies(num_ensembles);
	vector<vector2> pressures(num_ensembles);

	Energy energy; energy.E = E; energy.totE = totE; energy.sum_varE = sum_varE;
	cout << E << endl;

	vector<Count> counts(num_ensembles);
	vector<int> indices; vector<Jump> jumps(num_ensembles);
	Count temp_count;

	for (int i = 0; i < num_ensembles; ++i)
	{
		indices.push_back(i);
		string T_str; stringstream num;
		num << Temperatures[i]; num >> T_str;
		ofstream save("output_" + T_str + ".csv");
		states[i] = users;
		energies[i] = energy;
		pressures[i] = P;
		temp_count.accept_count = 0; temp_count.tot = 0, temp_count.tot_steps = 0;
		counts[i] = temp_count;
		jumps[i].dX = 0.1; jumps[i].dY = 0.1; jumps[i].dA = 5; 
		threads.push_back(std::thread(update, std::ref(states[i]), Temperatures[i], std::ref(energies[i]),  std::ref(pressures[i]), std::ref(counts[i]), std::ref(jumps[i]), 0));
	}
	for (thread &t : threads)
	{
	    if (t.joinable()) t.join();
	}
	double dbeta, dE;
	int index_i, index_j, temp_index;
	int swap_count = 0; int tot_count = 0;
	ofstream save0("check.dat");

	for (int k = 0; k < MAX_SWAPS; ++k)
	{
		for (int i = k%2; i+1 < num_ensembles; i += 2)
		{
			tot_count ++;
			index_i = indices[i]; index_j = indices[i+1];
			dbeta = 1/Temperatures[i] - 1/Temperatures[i+1];   /// temperatures never get swapped
			dE = energies[i].E - energies[i+1].E;
			if (SwapProb(dbeta, dE))
			{
				swap_count ++;
				temp_index = index_i; indices[i] = index_j; indices[i+1] = temp_index;  // swap configuration indices

				E = energies[i].E; energies[i].E = energies[i+1].E; energies[i+1].E = E;  // swap current energies

				pv = pressures[i].val2; pressures[i].val2 = pressures[i+1].val2; pressures[i+1].val2 = pv;
			}
		}
		if (k % 500 == 0) save0 << k << '\t' << 1.0*swap_count/tot_count << endl;

		int ensemble_index;
		for (int num = 0; num < num_ensembles; ++num)
		{
			ensemble_index = indices[num];
			threads[num] = std::thread(update, std::ref(states[ensemble_index]), Temperatures[num],  std::ref(energies[num]), std::ref(pressures[num]), std::ref(counts[num]), std::ref(jumps[num]), k);
		} // update each ensemble, keep temperature unchange, only change the configurations i.e., locations, energies, magnetizations, pressures becauce those calculations don't depend on T
		for (std::thread &t : threads)
		{
		    if (t.joinable()) t.join();
		}
	}
	cout << "Completed!" << std::endl;
};
