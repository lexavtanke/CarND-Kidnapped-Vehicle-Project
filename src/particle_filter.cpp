/**
 * particle_filter.cpp
 *
 * Created on: Dec 12, 2016
 * Author: Tiffany Huang
 */

#include "particle_filter.h"

#include <math.h>
#include <algorithm>
#include <iostream>
#include <iterator>
#include <numeric>
#include <random>
#include <string>
#include <vector>

#include "helper_functions.h"

using std::string;
using std::vector;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
  /**
   * TODO: Set the number of particles. Initialize all particles to 
   *   first position (based on estimates of x, y, theta and their uncertainties
   *   from GPS) and all weights to 1. 
   * TODO: Add random Gaussian noise to each particle.
   * NOTE: Consult particle_filter.h for more information about this method 
   *   (and others in this file).
   */
  // try to use as min number of particle as can 
  // but with 100 particles precision is 2 times better
  num_particles = 10;  // TODO: Set the number of particles
  
  std::default_random_engine gen;

  // normal distribution for x,y,theta

  std::normal_distribution<double> dist_x(x, std[0]);
  std::normal_distribution<double> dist_y(y, std[1]);
  std::normal_distribution<double> dist_theta(theta, std[2]);
  
  for (int i = 0; i < num_particles; ++i)
  {
    Particle sample;
    sample.id = i;
    sample.x = dist_x(gen);
    sample.y = dist_y(gen);
    sample.theta = dist_theta(gen);
    sample.weight = 1.0;
    particles.push_back(sample);
    weights.push_back(sample.weight);
  
  }
  is_initialized = true;
  // std::cout << std::endl << "num of particles is " << particles.size() << std::endl;

}

void ParticleFilter::prediction(double delta_t, double std_pos[], 
                                double velocity, double yaw_rate) {
  /**
   * TODO: Add measurements to each particle and add random Gaussian noise.
   * NOTE: When adding noise you may find std::normal_distribution 
   *   and std::default_random_engine useful.
   *  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
   *  http://www.cplusplus.com/reference/random/default_random_engine/
   */
  std::default_random_engine gen;

  for (int i = 0; i < particles.size(); ++i)
  {
      // define normal distributions for sensor noise
    std::normal_distribution<double> N_theta(0, std_pos[2]);
    std::normal_distribution<double> N_y(0, std_pos[1]);
    std::normal_distribution<double> N_x(0, std_pos[0]);

    // calculate new state
    // for working after the first circle
    if (fabs(yaw_rate) < 0.00001) {  
      particles[i].x += velocity * delta_t * cos(particles[i].theta);
      particles[i].y += velocity * delta_t * sin(particles[i].theta);
    } 
    else {
      particles[i].x += velocity / yaw_rate * (sin(particles[i].theta + yaw_rate*delta_t) - sin(particles[i].theta));
      particles[i].y += velocity / yaw_rate * (cos(particles[i].theta) - cos(particles[i].theta + yaw_rate*delta_t));
      particles[i].theta += yaw_rate * delta_t;
    }

    // add noise
    particles[i].x += N_x(gen);
    particles[i].y += N_y(gen);
    particles[i].theta += N_theta(gen);
  }

}

void ParticleFilter::dataAssociation(vector<LandmarkObs> predicted, 
                                     vector<LandmarkObs>& observations) {
  /**
   * TODO: Find the predicted measurement that is closest to each 
   *   observed measurement and assign the observed measurement to this 
   *   particular landmark.
   * NOTE: this method will NOT be called by the grading code. But you will 
   *   probably find it useful to implement this method and use it as a helper 
   *   during the updateWeights phase.
   */
  // for each observed measurement 
  for (int i = 0; i < observations.size(); ++i)
  {

    // set minimum distanse to big value
    double min_dist = std::numeric_limits<double>::max();
    // for each predicted measurement 
    for (int j = 0; j < predicted.size(); ++j)
    {
      // compute distanse between predicted and observed
      double distance = dist(predicted[j].x, predicted[j].y, observations[i].x, observations[i].y);
      // find the neares landmark
      if (distance < min_dist)
      {
        min_dist = distance;
        observations[i].id = predicted[j].id;
      }
    }
  }
}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
                                   const vector<LandmarkObs> &observations, 
                                   const Map &map_landmarks) {
  /**
   * TODO: Update the weights of each particle using a mult-variate Gaussian 
   *   distribution. You can read more about this distribution here: 
   *   https://en.wikipedia.org/wiki/Multivariate_normal_distribution
   * NOTE: The observations are given in the VEHICLE'S coordinate system. 
   *   Your particles are located according to the MAP'S coordinate system. 
   *   You will need to transform between the two systems. Keep in mind that
   *   this transformation requires both rotation AND translation (but no scaling).
   *   The following is a good resource for the theory:
   *   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
   *   and the following is a good resource for the actual equation to implement
   *   (look at equation 3.33) http://planning.cs.uiuc.edu/node99.html
   */
  //for each particle find all possible observed landmarks by sensor range 
  // and push back them to the prediction vector
  for (int i = 0; i < num_particles; ++i)
  {
    // coordinates of the particle
    double p_x = particles[i].x;
    double p_y = particles[i].y;
    double p_theta = particles[i].theta;
    // vector for possible observations 
    vector<LandmarkObs> predictions;
    // for each landmark in map
    for (int j = 0; j < map_landmarks.landmark_list.size(); ++j)
    {
      // get coords of the landmark
      float lm_x = map_landmarks.landmark_list[j].x_f;
      float lm_y = map_landmarks.landmark_list[j].y_f;
      int lm_id = map_landmarks.landmark_list[j].id_i;

      // if distanse from particle to landmark < sensor range add landmark to predictions
      if (dist(p_x, p_y, lm_x, lm_y) <= sensor_range)
      {
        predictions.push_back(LandmarkObs{lm_id, lm_x, lm_y});
      }
    }

      // tranform observations from vehicle coords to map 
    vector<LandmarkObs> trans_observ;
    for (int j = 0; j < observations.size(); ++j)
    {
      double trans_x = p_x + cos(p_theta) * observations[j].x - sin(p_theta) * observations[j].y;
      double trans_y = p_y + sin(p_theta) * observations[j].x + cos(p_theta) * observations[j].y;
      trans_observ.push_back(LandmarkObs{observations[j].id, trans_x, trans_y});
    } 

    // associate transformed observations with predictions
    // add landmarks id to observations
    dataAssociation(predictions, trans_observ);

    // updates weights
    // double new_weight = 1.0;
    particles[i].weight = 1.0;

    // std::cout << std::endl << "particle "<< i << std::endl;

    for (int j = 0; j < trans_observ.size(); ++j)
    {
      double obs_x = trans_observ[j].x;
      double obs_y = trans_observ[j].y;
      double pred_x;
      double pred_y;

      // find prediction of observation by id
      for (int n = 0; n < predictions.size(); n++)
      {
        if (predictions[n].id == trans_observ[j].id)
        {
          pred_x = predictions[n].x;
          pred_y = predictions[n].y; 
        }
      }
      // calculate weight for the observation
      double std_x = std_landmark[0];
      double std_y = std_landmark[1];
      double exponent = exp( -( pow(pred_x-obs_x,2)/(2*pow(std_x, 2)) + (pow(pred_y-obs_y,2)/(2*pow(std_y, 2))) ) );
      // std::cout << "exponent is " << exponent << std::endl;
      double normalizer = 2 * M_PI * std_x * std_y;
      // std::cout << "normalizer is " << normalizer << std::endl;
      double new_weight =  (1 / normalizer) * exponent;
      
      // std::cout << "current weight is " << particles[i].weight << " while new " << new_weight << std::endl;
      particles[i].weight *= new_weight;
    }


  }


}

void ParticleFilter::resample() {
  /**
   * TODO: Resample particles with replacement with probability proportional 
   *   to their weight. 
   * NOTE: You may find std::discrete_distribution helpful here.
   *   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
   */
  std::default_random_engine gen;

  vector<Particle> new_particles;

  // get all of the current weights
  vector<double> weights;
  for (int i = 0; i < num_particles; i++) {
    weights.push_back(particles[i].weight);
  }

  // generate random starting index for resampling wheel
  std::uniform_int_distribution<int> uniintdist(0, num_particles-1);
  auto index = uniintdist(gen);

  // get max weight
  double max_weight = *max_element(weights.begin(), weights.end());

  // uniform random distribution [0.0, max_weight)
  std::uniform_real_distribution<double> unirealdist(0.0, max_weight);

  double beta = 0.0;

  // spin the resample wheel!
  for (int i = 0; i < num_particles; i++) {
    beta += unirealdist(gen) * 2.0;
    while (beta > weights[index]) {
      beta -= weights[index];
      index = (index + 1) % num_particles;
    }
    new_particles.push_back(particles[index]);
  }

  particles = new_particles;
}

void ParticleFilter::SetAssociations(Particle& particle, 
                                     const vector<int>& associations, 
                                     const vector<double>& sense_x, 
                                     const vector<double>& sense_y) {
  // particle: the particle to which assign each listed association, 
  //   and association's (x,y) world coordinates mapping
  // associations: The landmark id that goes along with each listed association
  // sense_x: the associations x mapping already converted to world coordinates
  // sense_y: the associations y mapping already converted to world coordinates
  particle.associations= associations;
  particle.sense_x = sense_x;
  particle.sense_y = sense_y;
}

string ParticleFilter::getAssociations(Particle best) {
  vector<int> v = best.associations;
  std::stringstream ss;
  copy(v.begin(), v.end(), std::ostream_iterator<int>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length()-1);  // get rid of the trailing space
  return s;
}

string ParticleFilter::getSenseCoord(Particle best, string coord) {
  vector<double> v;

  if (coord == "X") {
    v = best.sense_x;
  } else {
    v = best.sense_y;
  }

  std::stringstream ss;
  copy(v.begin(), v.end(), std::ostream_iterator<float>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length()-1);  // get rid of the trailing space
  return s;
}