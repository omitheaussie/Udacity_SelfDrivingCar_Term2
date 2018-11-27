/*
* particle_filter.cpp
*
*  Created on: Dec 12, 2016
*      Author: Tiffany Huang
*/
//#include <random>
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

void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of 
	//   x, y, theta and their uncertainties from GPS) and all weights to 1. 
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).
	//default_random_engine gen;
	num_particles = 20;

	// This line creates a normal (Gaussian) distribution for x
	normal_distribution<double> dist_x(x, std[0]);

	// TODO: Create normal distributions for y and theta
	normal_distribution<double> dist_y(y, std[1]);
	normal_distribution<double> dist_theta(theta, std[2]);

	// Resize particles and weights variables
	weights.resize(num_particles);
	particles.resize(num_particles);

	// Initialize particles
	for (int i = 0; i < num_particles; i++) {
		particles[i].id = i;
		particles[i].x = dist_x(gen);
		particles[i].y = dist_y(gen);
		particles[i].theta = dist_theta(gen);
		particles[i].weight = 1.0;
		weights.push_back(1.0);
	}

	is_initialized = true;

}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/

	// Engine for later generation of particles
	

	// This line creates a normal (Gaussian) distribution for x
	normal_distribution<double> dist_x(0, std_pos[0]);

	// TODO: Create normal distributions for y and theta
	normal_distribution<double> dist_y(0, std_pos[1]);
	normal_distribution<double> dist_theta(0, std_pos[2]);

	// Different equations based on if yaw_rate is zero or not
	for (int i = 0; i < num_particles; i++) {

		if (fabs(yaw_rate) <= 0.00001) {
			particles[i].x += velocity * delta_t * cos(particles[i].theta) + dist_x(gen);
			particles[i].y += velocity * delta_t * sin(particles[i].theta) + dist_y(gen);
			// Theta remains unchanged
		}
		else {
			particles[i].x += (velocity / yaw_rate) * (sin(particles[i].theta + (yaw_rate * delta_t)) - sin(particles[i].theta)) + dist_x(gen);
			particles[i].y += (velocity / yaw_rate) * (cos(particles[i].theta) - cos(particles[i].theta + (yaw_rate * delta_t))) + dist_y(gen);
			particles[i].theta += yaw_rate * delta_t + dist_theta(gen);
		}
	}

}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.

	// For each in range landmark for each observation assign ID of landmark closest to the 
	// observation to observation ID.

	for (unsigned int i = 0; i < observations.size(); i++) { 
															
		double minDistance = numeric_limits<double>::max();

		// For each predition.
		for (int j = 0; j < predicted.size(); j++) { 

			double xDist = observations[i].x - predicted[j].x;
			double yDist = observations[i].y - predicted[j].y;

			double distance = sqrt(xDist * xDist + yDist * yDist);

			if (distance < minDistance || (minDistance == -1.0)) {
				minDistance = distance;
				observations[i].id = predicted[j].id;
			}
		}
	}

}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[],
	const std::vector<LandmarkObs> &observations, const Map &map_landmarks) {
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

	// Common terms in the multivariate Gaussian
	double common; 
	double const_common = 1. / (2.*M_PI*std_landmark[0] * std_landmark[1]);
	double x_denom = 2 * std_landmark[0] * std_landmark[0];
	double y_denom = 2 * std_landmark[1] * std_landmark[1];

	weights.clear();
	// Iterate for each particle for each landmark
	for (int i = 0; i < num_particles; i++) {

		particles[i].weight = 1.0; // Reinitialise weight to 1.0
		std::vector<LandmarkObs> trans_Obs;
		std::vector<LandmarkObs> pred_meas; // will only store landmarks that are within sensor range
		double cur_p_x = particles[i].x;
		double cur_p_y = particles[i].y;
		double cur_p_theta = particles[i].theta;

		// Landmarks which are within sensor range
		for (int k = 0; k < map_landmarks.landmark_list.size(); k++) {
			double cur_lm_x = map_landmarks.landmark_list[k].x_f;
			double cur_lm_y = map_landmarks.landmark_list[k].y_f;
			int id = map_landmarks.landmark_list[k].id_i;
			double d = dist(cur_p_x, cur_p_y, cur_lm_x, cur_lm_y);
			if (d < sensor_range) {
				LandmarkObs temp;
				temp = { id,cur_lm_x,cur_lm_y };
				pred_meas.push_back(temp);
			}
		}

		//Transform observations from vehicle to map coordinates
		for (int k = 0; k < observations.size(); k++) {

			double x;
			double y;
			int id;
			// Geometric Transformation
			x = observations[k].x*cos(cur_p_theta) - observations[k].y*sin(cur_p_theta) + cur_p_x;
			y = observations[k].x*sin(cur_p_theta) + observations[k].y*cos(cur_p_theta) + cur_p_y,
				id = observations[k].id;
			LandmarkObs temp;
			temp.x = x;
			temp.y = y;
			temp.id = id;
			trans_Obs.push_back(temp);

		}

		// Association
		dataAssociation(pred_meas, trans_Obs);
		//particles[i].weight = 1.0; // Reinitialise weight to 1.0
		std::vector<double> s_x;
		std::vector<double> s_y;
		std::vector<int> assoc;

		for (int n = 0; n < trans_Obs.size(); n++) {

			double mu_x = map_landmarks.landmark_list[trans_Obs[n].id - 1].x_f;
			double mu_y = map_landmarks.landmark_list[trans_Obs[n].id - 1].y_f;

			common = const_common * exp(-(((trans_Obs[n].x - mu_x)*(trans_Obs[n].x - mu_x)) / (2.*std_landmark[0] * std_landmark[0]) + ((trans_Obs[n].y - mu_y)*(trans_Obs[n].y - mu_y)) / (2.*std_landmark[1] * std_landmark[1])));
			particles[i].weight *= common;

			assoc.push_back(trans_Obs[n].id);
			s_x.push_back(trans_Obs[n].x);
			s_y.push_back(trans_Obs[n].y);
		}

		particles[i] = SetAssociations(particles[i], assoc, s_x, s_y);
		//cout << "ID" << i << "weight" << particles[i].weight << endl;
		weights.push_back(particles[i].weight);


	}
}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
	// Vector for new particles

	vector<Particle> resample_particles;
	discrete_distribution<int> distribution(weights.begin(), weights.end());
	for (int i = 0; i < num_particles; i++) {
		resample_particles.push_back(particles[distribution(gen)]);
	}

	particles = resample_particles;
}

Particle ParticleFilter::SetAssociations(Particle particle, const std::vector<int> associations,
	const std::vector<double> sense_x, const std::vector<double> sense_y)
{
	//particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
	// associations: The landmark id that goes along with each listed association
	// sense_x: the associations x mapping already converted to world coordinates
	// sense_y: the associations y mapping already converted to world coordinates

	particle.associations.clear();
	particle.sense_x.clear();
	particle.sense_y.clear();

	particle.associations = associations;
	particle.sense_x = sense_x;
	particle.sense_y = sense_y;

	return particle;
}

string ParticleFilter::getAssociations(Particle best)
{
	vector<int> v = best.associations;
	stringstream ss;
	copy(v.begin(), v.end(), ostream_iterator<int>(ss, " "));
	string s = ss.str();
	s = s.substr(0, s.length() - 1);  // get rid of the trailing space
	return s;
}
string ParticleFilter::getSenseX(Particle best)
{
	vector<double> v = best.sense_x;
	stringstream ss;
	copy(v.begin(), v.end(), ostream_iterator<float>(ss, " "));
	string s = ss.str();
	s = s.substr(0, s.length() - 1);  // get rid of the trailing space
	return s;
}
string ParticleFilter::getSenseY(Particle best)
{
	vector<double> v = best.sense_y;
	stringstream ss;
	copy(v.begin(), v.end(), ostream_iterator<float>(ss, " "));
	string s = ss.str();
	s = s.substr(0, s.length() - 1);  // get rid of the trailing space
	return s;
}
