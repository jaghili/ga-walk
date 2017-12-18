#include <iostream>
#include <algorithm>
#include <random>
#include <Eigen/Core>
#include <SFML/Graphics.hpp>
#include <SFML/Window.hpp>

const uint N = 100; // number of line nodes
const uint N_pop = 15; // number of species

typedef Eigen::Matrix<double, N, 2> Mat;

// initial position
const double xi = 10.;
const double yi = 10.;
//const Eigen::Vector2f vi(xi,yi);
// target position
const double xf = 700.;
const double yf = 700.;
//const Eigen::Vector2f vf(xf,yf);
const double L = std::sqrt(std::pow(xf-xi,2) + std::pow(yf-yi,2));

const double stepx = L/(N-1);
const double stepy = L/(N-1);
//const Eigen::Vector2f step(stepx, stepy);

const double tol = 0.1;

std::random_device rd;
std::mt19937 gen(rd());
std::uniform_real_distribution<double> dist_ur(-1.0,1.0);
std::uniform_int_distribution<> dist_ui(0,N-1);

// specie
struct steps {
  Mat p;
  double fitness = 1e300;
  
  //constructor
  steps() {
    // create steps
    p(0,0) = 0.;
    p(0,1) = 0.;
    for (int i = 1; i < N; i++) {            
      p(i,0) = dist_ur(gen) * stepx;
      p(i,1) = dist_ur(gen) * stepy;
    }  
  }
  
};

// fitness related functions
void computefit(steps & c) {
  c.fitness = std::sqrt(std::pow(xf-xi - c.p.col(0).sum(),2) + std::pow(yf- yi-c.p.col(1).sum(),2));
  //c.fitness = L - std::sqrt(std::pow((xi + c.p.col(0).sum()),2) + std::pow((yi + c.p.col(1).sum()),2));
}
bool orderbyfitness(steps & a, steps & b) { return (a.fitness < b.fitness); }

// mutation operator
void mutate(std::vector<steps> & pop) {
  // copy
  std::vector<steps> parent = pop; // copy
  // blank
  for (auto & offspring : pop) { offspring.p = Eigen::MatrixXd::Zero(N,2); }
  
  // Fill with offsprings
  // take the two best from pp
  pop[0].p = parent[0].p;
  // mean
  pop[1].p = 0.5 * parent[0].p + 0.5 * parent[1].p;

  // weighted sum
  pop[2].p = 0.1 * parent[4].p + 0.9 * parent[1].p;

  // crossover between 0 and 1
  pop[3].p = parent[0].p;
  pop[3].p.bottomLeftCorner(std::floor((double)N/2),2) = parent[1].p.bottomLeftCorner(std::floor((double)N/2),2);

  // random pointwise mutation
  pop[4].p = parent[0].p;
  pop[4].p(dist_ui(gen),0) = dist_ur(gen) * stepx;
  pop[4].p(dist_ui(gen),1) = dist_ur(gen) * stepy;
}

int main () {
  // Init population
  std::vector<steps> pop;
  pop.reserve(N_pop);
  for (int i = 0; i < N_pop; i++) {
    steps p;
    pop.push_back(p);
  }

  // Plot : Render window
  double global_fit = 1e300;    
  
  sf::RenderWindow window(sf::VideoMode(800,800), "SFML Plot");
  window.setFramerateLimit(30); // FPS limiter

  sf::CircleShape target(20);
  target.setPosition(sf::Vector2f(xf,yf));
  target.setFillColor(sf::Color::Red);

  uint iter = 0;
  while (global_fit > tol) {
    window.clear(sf::Color(0,0,0));

    // draw target
    window.draw(target);
    
    global_fit = 1e300;
    for (auto & s : pop) {
      // compute fitness scores
      computefit(s);
      global_fit = std::min(global_fit,s.fitness);

      // draw
      double _x = xi;
      double _y = yi;
      for (int i = 0; i < N; i++) {
	sf::Vertex line[] = {
	  sf::Vertex(sf::Vector2f( _x, _y), sf::Color::Red),
	  sf::Vertex(sf::Vector2f( _x + s.p(i,0), _y + s.p(i,1)), sf::Color::Red)
	};
	window.draw(line, 2, sf::Lines);
	_x = _x + s.p(i,0);
	_y = _y + s.p(i,1);
      }
    }
    window.display();

    // mutation
    if (global_fit > tol) {
      for (int i = 0; i < N_pop; i++) {	std::cout << "fit" << i << "=" << pop[i].fitness << " "; }
      std::cout << "\t global fit=" << global_fit << " -> mutating..." << std::endl;
      
      std::sort(pop.begin(), pop.end(), orderbyfitness); // reorder
      mutate(pop); // mutate
      iter++;
    }
    else {
      std::cout << "ACHIEVED! iters=" << iter<< " global_fit=" << global_fit << std::endl;
    }
  }
  // exit
  window.close();

  return 0;
}
