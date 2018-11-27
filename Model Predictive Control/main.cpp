#include <math.h>
#include <uWS/uWS.h>
#include <chrono>
#include <iostream>
#include <thread>
#include <vector>
#include "Eigen-3.3/Eigen/Core"
#include "Eigen-3.3/Eigen/QR"
#include "MPC.h"
#include "json.hpp"
#define DEBUG 1
// for convenience
using json = nlohmann::json;

// For converting back and forth between radians and degrees.
constexpr double pi() { return M_PI; }
double deg2rad(double x) { return x * pi() / 180; }
double rad2deg(double x) { return x * 180 / pi(); }

// Checks if the SocketIO event has JSON data.
// If there is data the JSON object in string format will be returned,
// else the empty string "" will be returned.
string hasData(string s) {
  auto found_null = s.find("null");
  auto b1 = s.find_first_of("[");
  auto b2 = s.rfind("}]");
  if (found_null != string::npos) {
    return "";
  } else if (b1 != string::npos && b2 != string::npos) {
    return s.substr(b1, b2 - b1 + 2);
  }
  return "";
}

// Evaluate a polynomial.
double polyeval(Eigen::VectorXd coeffs, double x) {
  double result = 0.0;
  for (int i = 0; i < coeffs.size(); i++) {
    result += coeffs[i] * pow(x, i);
  }
  return result;
}

// Fit a polynomial.
// Adapted from
// https://github.com/JuliaMath/Polynomials.jl/blob/master/src/Polynomials.jl#L676-L716
Eigen::VectorXd polyfit(Eigen::VectorXd xvals, Eigen::VectorXd yvals,
                        int order) {
  assert(xvals.size() == yvals.size());
  assert(order >= 1 && order <= xvals.size() - 1);
  Eigen::MatrixXd A(xvals.size(), order + 1);

  for (int i = 0; i < xvals.size(); i++) {
    A(i, 0) = 1.0;
  }

  for (int j = 0; j < xvals.size(); j++) {
    for (int i = 0; i < order; i++) {
      A(j, i + 1) = A(j, i) * xvals(j);
    }
  }

  auto Q = A.householderQr();
  auto result = Q.solve(yvals);
  return result;
}

int main() {
  uWS::Hub h;

  // MPC is initialized here!
  MPC mpc;

  h.onMessage([&mpc](uWS::WebSocket<uWS::SERVER> ws, char *data, size_t length,
                     uWS::OpCode opCode) {
    // "42" at the start of the message means there's a websocket message event.
    // The 4 signifies a websocket message
    // The 2 signifies a websocket event
    string sdata = string(data).substr(0, length);
    cout << sdata << endl;
    if (sdata.size() > 2 && sdata[0] == '4' && sdata[1] == '2') {
      string s = hasData(sdata);
      if (s != "") {
        auto j = json::parse(s);
        string event = j[0].get<string>();
        if (event == "telemetry") {
          // j[1] is the data JSON object
		  //"ptsx":[-32.16173,-43.49173,-61.09,-78.29172,-93.05002,-107.7717]
		  //"ptsy":[113.361,105.941,92.88499,78.73102,65.34102,50.57938]
		  //"psi_unity":4.12033,"psi":3.733651,"x":-40.62,"y":108.73
		  //"steering_angle":0
		  //"throttle":0
		  //"speed":0.4380091
		  
          vector<double> ptsx = j[1]["ptsx"];
          vector<double> ptsy = j[1]["ptsy"];
          double px = j[1]["x"];
          double py = j[1]["y"];
          double psi = j[1]["psi"];
          double v = j[1]["speed"];
          double delta= j[1]["steering_angle"];
          double throttle = j[1]["throttle"];

          /*
          * TODO: Calculate steering angle and throttle using MPC.
          *
          * Both are in between [-1, 1].
          *
          */
		#if(DEBUG)
			std::cout << "Calculating steering angle" << std::endl;
		#endif
          // Transforms waypoints global coordinates to the vehicle coordinates.
          auto ptx = Eigen::VectorXd(ptsx.size());
		  auto pty = Eigen::VectorXd(ptsy.size()); 
          for (int i = 0; i < ptsx.size(); i++) {
			// Subtract px,py from ptsx,ptsy as during transformations vehicle position px,py is 
			// the origin and then transform the rest of the ptsx,ptsy using the psi to 
			// rotate the coordinates in veh coordinates
            double tempx = ptsx[i] - px;
            double tempy = ptsy[i] - py;
            double psi_rad = psi;
            ptx[i] = tempx*cos(0 - psi_rad) - tempy*sin(0 - psi_rad);
            pty[i] = tempx*sin(0 - psi_rad) + tempy*cos(0 - psi_rad);
          }

          // Fit polynomial to the points - 3rd order.
		  // This is basically drawing a fitted curve between the transformed coordinates
		  int poly_order = 3;
          auto coeffs = polyfit(ptx, pty, poly_order); // coeffs provides 3rd order function 
													   // with coeffs of each order in "coeffs" variable
		  
		  // By evaluating the above function at 0,0 (px,py - origin) we find the cte
		  // So cte can actually be the first coeff
		  double cte = -polyeval(coeffs, 0.);    
		  // Derivative of the function gives its slope in radians which is error_psi = epsi
          double epsi = -atan(coeffs[1]);
		  
          // Actuator delay in milliseconds.
          const int actuatorDelay =  100;

          // Actuator delay in seconds.
          const double delay = actuatorDelay / 1000.0;
		
          // Initial state.
          const double x0 = 0;
          const double y0 = 0;
          const double psi0 = 0;
          const double cte0 = coeffs[0];
          const double epsi0 = -atan(coeffs[1]);

          // State after delay.
          double x_delay = x0 + ( v * cos(psi0) * delay );
          double y_delay = y0 + ( v * sin(psi0) * delay );
          double psi_delay = psi0 - ( v * delta * delay / mpc.Lf );
          double v_delay = v + throttle * delay;
          double cte_delay = cte0 + ( v * sin(epsi0) * delay );
          double epsi_delay = epsi0 - ( v * atan(coeffs[1]) * delay / mpc.Lf );

          // Define the state vector.
          Eigen::VectorXd state(6);
          state << x_delay, y_delay, psi_delay, v_delay, cte_delay, epsi_delay;

          // Find the MPC solution.
          auto vars = mpc.Solve(state, coeffs);

          double steer_value = vars[0]/deg2rad(25);
          double throttle_value = vars[1];

          json msgJson;
          // NOTE: Remember to divide by deg2rad(25) before you send the steering value back.
          // Otherwise the values will be in between [-deg2rad(25), deg2rad(25] instead of [-1, 1].
          msgJson["steering_angle"] = steer_value;
          msgJson["throttle"] = throttle_value;

          //Display the MPC predicted trajectory
          vector<double> mpc_x_vals;
          vector<double> mpc_y_vals;

          for ( int i = 2; i < vars.size(); i++ ) {
            if ( i % 2 == 0 ) {
              mpc_x_vals.push_back( vars[i] );
            } else {
              mpc_y_vals.push_back( vars[i] );
            }
          }

          //.. add (x,y) points to list here, points are in reference to the vehicle's coordinate system
          // the points in the simulator are connected by a Green line

          msgJson["mpc_x"] = mpc_x_vals;
          msgJson["mpc_y"] = mpc_y_vals;

          //Display the waypoints/reference line
          vector<double> next_x_vals;
          vector<double> next_y_vals;

          int n_points = 5; //10
          int step_size = 8;
          for (int i = 0; i < n_points; i++) {
            double xstep = (double)std::pow(1.5, i)*step_size;
            next_x_vals.push_back(xstep);
            next_y_vals.push_back(polyeval(coeffs, xstep));
          }

          //.. add (x,y) points to list here, points are in reference to the vehicle's coordinate system
          // the points in the simulator are connected by a Yellow line

          msgJson["next_x"] = next_x_vals;
          msgJson["next_y"] = next_y_vals;


          auto msg = "42[\"steer\"," + msgJson.dump() + "]";
          // std::cout << msg << std::endl;
          // Latency
          // The purpose is to mimic real driving conditions where
          // the car does actuate the commands instantly.
          //
          // Feel free to play around with this value but should be to drive
          // around the track with 100ms latency.
          //
          // NOTE: REMEMBER TO SET THIS TO 100 MILLISECONDS BEFORE
          // SUBMITTING.
          this_thread::sleep_for(chrono::milliseconds(actuatorDelay));
          ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
        }
      } else {
        // Manual driving
        std::string msg = "42[\"manual\",{}]";
        ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
      }
    }
  });

  // We don't need this since we're not using HTTP but if it's removed the
  // program
  // doesn't compile :-(
  h.onHttpRequest([](uWS::HttpResponse *res, uWS::HttpRequest req, char *data,
                     size_t, size_t) {
    const std::string s = "<h1>Hello world!</h1>";
    if (req.getUrl().valueLength == 1) {
      res->end(s.data(), s.length());
    } else {
      // i guess this should be done more gracefully?
      res->end(nullptr, 0);
    }
  });

  h.onConnection([&h](uWS::WebSocket<uWS::SERVER> ws, uWS::HttpRequest req) {
    std::cout << "Connected!!!" << std::endl;
  });

  h.onDisconnection([&h](uWS::WebSocket<uWS::SERVER> ws, int code,
                         char *message, size_t length) {
    ws.close();
    std::cout << "Disconnected" << std::endl;
  });

  int port = 4567;
  auto host = "127.0.0.1";
  if (h.listen(host,port)) {
    std::cout << "Listening to port " << port << std::endl;
  } else {
    std::cerr << "Failed to listen to port" << std::endl;
    return -1;
  }
  h.run();
}
