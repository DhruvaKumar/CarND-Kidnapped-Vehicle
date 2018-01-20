/*
 * particle_filter.cpp
 *
 *  Created on: Dec 12, 2016
 *      Author: Tiffany Huang
 */

#include <random>
#include <algorithm>
#include <iostream>
#include <numeric>
#include <math.h> 
#include <iostream>
#include <sstream>
#include <string>
#include <iterator>

#include "particle_filter.h"

using namespace std;

static const double EPSILON = 0.00001;

void ParticleFilter::init(double x, double y, double theta, double std[]) 
{
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of 
	//   x, y, theta and their uncertainties from GPS) and all weights to 1. 
	// Add random Gaussian noise to each particle.

	// initialize particles. sample from gaussian distribution centered around gps estimates
	std::normal_distribution<> x_dist{x, std[0]};
	std::normal_distribution<> y_dist{y, std[1]};
	std::normal_distribution<> theta_dist{theta, std[2]};

	particles.clear();
	for (int i = 0; i < num_particles; ++i)
	{
		particles.emplace_back(x_dist(gen_), y_dist(gen_), theta_dist(gen_), 1.0);
	}

	std::cout << num_particles << " particles initialized around gps " << x << "," << y << "," << theta << "\n";
	// // debug
	// for (const auto& p : particles)
	// 	std::cout << "p:" << p.x << "," << p.y << "," << p.theta << " w:" << p.weight << "\n";

	is_initialized = true;
}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) 
{
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/

	// apply odometery + gaussian noise to every particle
	std::normal_distribution<> x_dist{0, std_pos[0]};
	std::normal_distribution<> y_dist{0, std_pos[1]};
	std::normal_distribution<> theta_dist{0, std_pos[2]};

	// // debug
	// std::cout << "pred: v:" << velocity << " yaw_rate:" << yaw_rate << "\n"; 

	for (auto& p : particles)
	{
		if (fabs(yaw_rate) < EPSILON)
		{
			p.x += velocity*delta_t*cos(p.theta) + x_dist(gen_);
			p.y += velocity*delta_t*sin(p.theta) + y_dist(gen_);
			p.theta += theta_dist(gen_);
		}
		else
		{
			p.x += velocity/yaw_rate * (sin(p.theta + yaw_rate*delta_t) - sin(p.theta)) + x_dist(gen_);
			p.y += velocity/yaw_rate * (cos(p.theta) - cos(p.theta + yaw_rate*delta_t)) + y_dist(gen_);
			p.theta += yaw_rate*delta_t + theta_dist(gen_);
		}

		// // debug
		// std::cout << &p-&particles[0] << ":" << p.x << "," << p.y << "," << p.theta << "\n";
	}
}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.

}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
		const std::vector<LandmarkObs> &observations, const Map &map_landmarks) 
{
	// TODO: Update the weights of each particle using a mult-variate Gaussian distribution. You can read
	//   more about this distribution here: https://en.wikipedia.org/wiki/Multivariate_normal_distribution
	// NOTE: The observations are given in the VEHICLE'S coordinate system. Your particles are located
	//   according to the MAP'S coordinate system. You will need to transform between the two systems.
	//   Keep in mind that this transformation requires both rotation AND translation (but no scaling).
	//   The following is a good resource for the theory:
	//   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
	//   and the following is a good resource for the actual equation to implement (look at equation 
	//   3.33
	//   http://planning.cs.uiuc.edu/node99.html

	double sum_weights = 0.0;

	for (auto& p : particles)
	{
		double weight = 1.0;
		
		// [optional] association graphics
		// p.associations.clear();
		// p.sense_x.clear();
		// p.sense_y.clear();
		// p.associations.reserve(observations.size());
		// p.sense_x.reserve(observations.size());
		// p.sense_y.reserve(observations.size());
		

		for (const auto& z : observations)
		{
			// transform observation z from vehicle/particle's coordinate system to the global map's coordinate system
			auto zx_m = p.x + z.x*cos(p.theta) - z.y*sin(p.theta);
			auto zy_m = p.y + z.x*sin(p.theta) + z.y*cos(p.theta);  

			// associate obs z with nearest landmark
			// TODO: scope for optimization (use kd tree to reduce search to O(lgn))
			double shortest_dist = numeric_limits<double>::max();
			LandmarkObs nearest_landmark;
			for (const auto& l : map_landmarks.landmark_list)
			{
				auto dist_obs_landmk = dist(zx_m, zy_m, l.x_f, l.y_f);
				if (dist_obs_landmk < shortest_dist)
				{
					nearest_landmark.id = l.id_i;
					nearest_landmark.x = l.x_f;
					nearest_landmark.y = l.y_f;
					shortest_dist = dist_obs_landmk;
				}
			}

			// // debug
			// auto update_eq = 1.0/(2*M_PI*std_landmark[0]*std_landmark[1]) * 
			// 	exp(-(pow((zx_m - nearest_landmark.x), 2)/(2*pow(std_landmark[0],2)) +
			// 		pow((zy_m - nearest_landmark.y),2)/(2*pow(std_landmark[1],2))));

			// // debug
			// std::cout << "updatew: " << update_eq << " x:" << zx_m << "," << nearest_landmark.x << " y:" << zy_m << "," << nearest_landmark.y << " p:" << p.x << "," << p.y << "\n";

			// update weight: P(zk | xt, M)
			weight *= 1.0/(2*M_PI*std_landmark[0]*std_landmark[1]) * 
				exp(-(pow((zx_m - nearest_landmark.x), 2)/(2*pow(std_landmark[0],2)) +
					pow((zy_m - nearest_landmark.y),2)/(2*pow(std_landmark[1],2))));

			// [optional] set associations
			// p.associations.push_back(nearest_landmark.id);
			// p.sense_x.push_back(nearest_landmark.x);
			// p.sense_y.push_back(nearest_landmark.y);
		}

		// if no observation detected by particle, set low weight 
		if (weight == 1.0)
		{
			std::cout << "no observations detected. setting low weight\n";
			weight = 1.0e-20;
		}

		// retain prior weight info from previous time step
		p.weight *= weight;

		sum_weights += p.weight;
	}

	// // debug
	// std::cout << "before norm ";
	// for (const auto& p : particles)
	// 	std::cout << p.weight << " ";
	// std::cout << "\n";


	// normalize weights
	std::cout << "weights: ";
	for (auto& p : particles)
	{
		p.weight /= sum_weights;

		// debug
		std::cout << p.weight << " ";
	}
	std::cout << "\n";
}

void ParticleFilter::resample() 
{
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution

	std::vector<double> weights;
	weights.reserve(particles.size());
	for (const auto& p : particles)
		weights.push_back(p.weight);

	std::discrete_distribution<> d(weights.begin(), weights.end());
	
	auto particles_before = particles;

	for (int i = 0; i < num_particles; ++i)
		particles[i] = particles_before[d(gen_)];
	
	// // debug
	// std::cout << "weights after resampling: ";
	// for (const auto& p : particles)
	// 	std::cout << p.weight << " ";
	// std::cout << "\n";

}

Particle ParticleFilter::SetAssociations(Particle& particle, const std::vector<int>& associations, 
                                     const std::vector<double>& sense_x, const std::vector<double>& sense_y)
{
    //particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
    // associations: The landmark id that goes along with each listed association
    // sense_x: the associations x mapping already converted to world coordinates
    // sense_y: the associations y mapping already converted to world coordinates

    particle.associations= associations;
    particle.sense_x = sense_x;
    particle.sense_y = sense_y;
}

string ParticleFilter::getAssociations(Particle best)
{
	vector<int> v = best.associations;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<int>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseX(Particle best)
{
	vector<double> v = best.sense_x;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseY(Particle best)
{
	vector<double> v = best.sense_y;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
