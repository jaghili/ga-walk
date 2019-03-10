#include<iostream>
#include<algorithm>
#include<random>
#include<Eigen/Core>
#include<SFML/Graphics.hpp>
#include<SFML/Window.hpp>

const uint n_points = 300;
const uint i_crossover = 299; // must be < n_points !!!
const uint n_species = 50;

typedef Eigen::Matrix<double, n_points, 2> Mat;
typedef Eigen::Vector2d Point;

// Window size
const double wsizex = 1000;
const double wsizey = 800;

// Initial point
const double xi = 10.0;
const double yi = 10.0;

// Target point
const double xf = 400.0;
const double yf = 600.0;
const double rf = 10;
sf::CircleShape target(rf);

// Obstacle
const double xo = 400.0;
const double yo = 400.0;
const double ro = 200.0;
double fobs(double x, double y) {
  //return (std::sqrt(std::pow((x - xo),2) + std::pow((y - yo),2)) < ro);
  return (std::sin(50*x/wsizex)*std::sin(50*y/wsizey)>0.7);
}

void draw_obstacle(sf::RenderWindow &w) {
  for (double _x = 0; _x < wsizex; _x += 5) {
    for (double _y = 0; _y < wsizey; _y += 5) {     
      if (fobs(_x,_y) > 0.) {
	sf::CircleShape cercle(2);
	cercle.setPosition(sf::Vector2f(_x,_y));
	cercle.setFillColor(sf::Color::White);
	w.draw(cercle);
      }
    }
  }
}


// steps
const double L = std::sqrt((xf-xi)*(xf-xi) + (yf-yi)*(yf-yi));
const double stepx = 1.5*L/n_points;
const double stepy = 1.5*L/n_points;
const double tgv = 1e30;

// random generators
std::random_device rd;
std::mt19937 gen(rd());
std::uniform_real_distribution<double> dist_ur(-1.0,1.0); // real 
std::uniform_int_distribution<> dist_ui(0,n_points-1); // integers


struct path {
  Mat points;
  double fitness = 1e300;

  // ctor
  path () {
    points(0,0) = xi;
    points(0,1) = yi;
    for (uint i = 1; i < n_points; i++) {
      points(i,0) = dist_ur(gen)*stepx;
      points(i,1) = dist_ur(gen)*stepy;
    }    
  }

  // draw the path
  void plot(sf::RenderWindow & w) {
    // init point
    double _x = xi;
    double _y = yi;
    
    for (uint i = 0; i < n_points; i++) {
      // init
      double dx  = points(i,0);
      double dy  = points(i,1);
      
      // Draw red line
      sf::Vertex line[] =
	{
	 sf::Vertex(sf::Vector2f(_x, _y), sf::Color::Red),
	 sf::Vertex(sf::Vector2f(_x + dx, _y + dy), sf::Color::Red)
	};
      w.draw(line, 2, sf::Lines);

      /*
      sf::CircleShape cercle(5);
      cercle.setPosition(sf::Vector2f(_x,_y));
      cercle.setFillColor(sf::Color::Blue);
      w.draw(cercle);
      */
      
      // update
      _x += dx;
      _y += dy;
    }
  }
};

void computeFitness(path & c) {
  double _x = xi;
  double _y = yi;
  c.fitness = 0.;
  for (uint i = 0; i < n_points; i++) {
    _x += c.points(i,0);
    _y += c.points(i,1);
    c.fitness += tgv * fobs(_x,_y); // penalize in case of crossing an obstacle
  }
  c.fitness += std::sqrt( (_x-xf)*(_x-xf) + (_y - yf)*(_y - yf) );
}


void mutate (std::vector<path> & chemins) {
  // copy 
  std::vector<path> parents = chemins;

  // reinit
  for (auto & c : chemins) { c.points = Eigen::MatrixXd::Zero(n_points,2); };

  chemins[0].points = parents[0].points;

  chemins[1].points = 0.9*parents[0].points + 0.1*parents[1].points;

  chemins[2].points = 0.5*parents[1].points + 0.5*parents[2].points;

  // crossover between 0 and 1
  chemins[3].points = parents[0].points;
  chemins[3].points.bottomLeftCorner(n_points-i_crossover,2) = parents[1].points.bottomLeftCorner(n_points-i_crossover,2);

  // crossover between 0 and 2
  uint r_cross = dist_ui(gen);
  chemins[4].points = parents[0].points;
  chemins[4].points.bottomLeftCorner(n_points-r_cross,2) = parents[2].points.bottomLeftCorner(n_points-r_cross,2);


  chemins[5].points = 0.7*parents[0].points + 0.3*parents[2].points;

  chemins[6].points = 0.6*parents[1].points + 0.4*parents[2].points;

  chemins[7].points = parents[0].points;
  chemins[7].points(n_points-1,0) = dist_ur(gen) * stepx;
  chemins[7].points(n_points-1,1) = dist_ur(gen) * stepy;
  
  for (uint i = 8; i < n_species; i++) {
    chemins[i].points = parents[0].points;
    uint mut_point = dist_ui(gen);    
    chemins[i].points(mut_point,0) = dist_ur(gen) * stepx;
    chemins[i].points(mut_point,1) = dist_ur(gen) * stepy;
  }
}




bool reorder(path & a, path & b) {
  return (a.fitness < b.fitness);
}



int main () {
  // Window Init
  sf::RenderWindow window(sf::VideoMode(wsizex, wsizey), "SFML Plot");
  //window.setFramerateLimit(60); // FPS limiter
  target.setPosition(sf::Vector2f(xf,yf));
  target.setFillColor(sf::Color::Red);

  double gfit = 1e300;
  uint n_gen = 0;
  
  // init chemins
  std::cout << "init paths..." << std::endl;
  std::vector<path> chemins;
  chemins.reserve(n_species);
  for (uint i = 0; i < n_species; i++) {
    path c;
    computeFitness(c);
    gfit = std::min(c.fitness, gfit);
    chemins.push_back(c);
  }


  // Loop 
  while (gfit > rf) {
    n_gen++;

    // reorder
    std::sort(chemins.begin(), chemins.end(), reorder);

    // plot
    window.clear(sf::Color(0,0,0));
    draw_obstacle(window);
    window.draw(target);
    chemins[0].plot(window); // plot only the best 
    //for (auto & c : chemins) { c.plot(window); } // plot everybody
    window.display();    

    // mutate
    mutate(chemins);

    // recompute fitness
    gfit = 1e300;
    for (auto & c : chemins) {
      computeFitness(c);
      gfit = std::min(c.fitness, gfit);
    }
    std::cout << n_gen << "\t" << gfit << std::endl;
  }

  // converged, plotting
  std::sort(chemins.begin(), chemins.end(), reorder);
  path & cbest = chemins[0];
  std::cout << "Converged in " << n_gen << " iters! fit=" << cbest.fitness <<  " plotting..." << std::endl;
  
  while (true) {
    window.clear(sf::Color(0,0,0));
    draw_obstacle(window);
    window.draw(target);
    cbest.plot(window);
    window.display();    
  }
  
  return 0;
}
