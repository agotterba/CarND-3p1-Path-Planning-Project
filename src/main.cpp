#include <fstream>
#include <math.h>
#include <uWS/uWS.h>
#include <chrono>
#include <iostream>
#include <thread>
#include <vector>
#include "Eigen-3.3/Eigen/Core"
#include "Eigen-3.3/Eigen/QR"
#include "json.hpp"
#include "spline.h"
#include "Eigen-3.3/Eigen/LU"

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;

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
  auto b2 = s.find_first_of("}");
  if (found_null != string::npos) {
    return "";
  } else if (b1 != string::npos && b2 != string::npos) {
    return s.substr(b1, b2 - b1 + 2);
  }
  return "";
}

double distance(double x1, double y1, double x2, double y2)
{
	return sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1));
}
int ClosestWaypoint(double x, double y, const vector<double> &maps_x, const vector<double> &maps_y)
{

	double closestLen = 100000; //large number
	int closestWaypoint = 0;

	for(int i = 0; i < maps_x.size(); i++)
	{
		double map_x = maps_x[i];
		double map_y = maps_y[i];
		double dist = distance(x,y,map_x,map_y);
		if(dist < closestLen)
		{
			closestLen = dist;
			closestWaypoint = i;
		}

	}

	return closestWaypoint;

}

int NextWaypoint(double x, double y, double theta, const vector<double> &maps_x, const vector<double> &maps_y)
{

	int closestWaypoint = ClosestWaypoint(x,y,maps_x,maps_y);

	double map_x = maps_x[closestWaypoint];
	double map_y = maps_y[closestWaypoint];

	double heading = atan2((map_y-y),(map_x-x));

	double angle = fabs(theta-heading);
  angle = min(2*pi() - angle, angle);

  if(angle > pi()/4)
  {
    closestWaypoint++;
  if (closestWaypoint == maps_x.size())
  {
    closestWaypoint = 0;
  }
  }

  return closestWaypoint;
}

// Transform from Cartesian x,y coordinates to Frenet s,d coordinates
vector<double> getFrenet(double x, double y, double theta, const vector<double> &maps_x, const vector<double> &maps_y)
{
	int next_wp = NextWaypoint(x,y, theta, maps_x,maps_y);

	int prev_wp;
	prev_wp = next_wp-1;
	if(next_wp == 0)
	{
		prev_wp  = maps_x.size()-1;
	}

	double n_x = maps_x[next_wp]-maps_x[prev_wp];
	double n_y = maps_y[next_wp]-maps_y[prev_wp];
	double x_x = x - maps_x[prev_wp];
	double x_y = y - maps_y[prev_wp];

	// find the projection of x onto n
	double proj_norm = (x_x*n_x+x_y*n_y)/(n_x*n_x+n_y*n_y);
	double proj_x = proj_norm*n_x;
	double proj_y = proj_norm*n_y;

	double frenet_d = distance(x_x,x_y,proj_x,proj_y);

	//see if d value is positive or negative by comparing it to a center point

	double center_x = 1000-maps_x[prev_wp];
	double center_y = 2000-maps_y[prev_wp];
	double centerToPos = distance(center_x,center_y,x_x,x_y);
	double centerToRef = distance(center_x,center_y,proj_x,proj_y);

	if(centerToPos <= centerToRef)
	{
		frenet_d *= -1;
	}

	// calculate s value
	double frenet_s = 0;
	for(int i = 0; i < prev_wp; i++)
	{
		frenet_s += distance(maps_x[i],maps_y[i],maps_x[i+1],maps_y[i+1]);
	}

	frenet_s += distance(0,0,proj_x,proj_y);

	return {frenet_s,frenet_d};

}

// Transform from Frenet s,d coordinates to Cartesian x,y
vector<double> getXY(double s, double d, const vector<double> &maps_s, const vector<double> &maps_x, const vector<double> &maps_y)
{
	int prev_wp = -1;

	while(s > maps_s[prev_wp+1] && (prev_wp < (int)(maps_s.size()-1) ))
	{
		prev_wp++;
	}

	int wp2 = (prev_wp+1)%maps_x.size();

	double heading = atan2((maps_y[wp2]-maps_y[prev_wp]),(maps_x[wp2]-maps_x[prev_wp]));
	// the x,y,s along the segment
	double seg_s = (s-maps_s[prev_wp]);

	double seg_x = maps_x[prev_wp]+seg_s*cos(heading);
	double seg_y = maps_y[prev_wp]+seg_s*sin(heading);

	double perp_heading = heading-pi()/2;

	double x = seg_x + d*cos(perp_heading);
	double y = seg_y + d*sin(perp_heading);

	return {x,y};

}

vector<double> jmt(vector<double> start,vector<double> end, double T){


  vector <double> a(6);
  a[0] = start[0];
  a[1] = start[1];
  a[2] = 0.5 * start[2];
  VectorXd eq(3);
  eq[0] = end[0] - (start[0] + start[1] * T + 0.5 *start[2] * T * T);
  eq[1] = end[1] - (start[1] + start[2] * T);
  eq[2] = end[2] - start[2];
  
  MatrixXd tmat(3,3);
  tmat <<            pow(T,3),   pow(T,4),    pow(T,5),
                   3*pow(T,2), 4*pow(T,3),  5*pow(T,4),
                   6*T,       12*pow(T,2), 20*pow(T,3);
    
  VectorXd x(3);
  //x = tmat.colPivHouseholderQr().solve(eq);
  MatrixXd tmati = tmat.inverse();
  x = tmati*eq;
  
  a[3] = x[0];
  a[4] = x[1];
  a[5] = x[2];
  
  return a;
}

vector<double> compute_lane_change(int old_lane, int new_lane, int lane_width,int num_lane_change_points){
  vector<double> start_d {(old_lane + 0.5) * lane_width,0.0,0.0};
  vector<double> end_d   {(new_lane + 0.5) * lane_width,0.0,0.0};
  vector<double> coeffs = jmt(start_d,end_d,num_lane_change_points);
  vector<double> lane_change_points;
  for (int i = 0; i < num_lane_change_points; i++){
    double dval = 0.0;
    double t = 1.0;
    for (int j = 0; j < coeffs.size(); j++){
      dval += t * coeffs[j];
      t *= i;
    }
    lane_change_points.push_back(dval);
  }
  return lane_change_points;
}

int main() {
  uWS::Hub h;

  // Load up map values for waypoint's x,y,s and d normalized normal vectors
  vector<double> map_waypoints_x;
  vector<double> map_waypoints_y;
  vector<double> map_waypoints_s;
  vector<double> map_waypoints_dx;
  vector<double> map_waypoints_dy;

  // Waypoint map to read from
  string map_file_ = "../data/highway_map.csv";
  // The max s value before wrapping around the track back to 0
  double max_s = 6945.554;

  ifstream in_map_(map_file_.c_str(), ifstream::in);

  string line;
  while (getline(in_map_, line)) {
  	istringstream iss(line);
  	double x;
  	double y;
  	float s;
  	float d_x;
  	float d_y;
  	iss >> x;
  	iss >> y;
  	iss >> s;
  	iss >> d_x;
  	iss >> d_y;
  	map_waypoints_x.push_back(x);
  	map_waypoints_y.push_back(y);
  	map_waypoints_s.push_back(s);
  	map_waypoints_dx.push_back(d_x);
  	map_waypoints_dy.push_back(d_y);
  }

  double ref_vel  = 0.0; //velocity for reference point.  declare ref_vel here so that it's persistent across calls to onMessage; start with car stopped
  int old_lane    = 1;   //current lane, or lane we're leaving.  lanes are 0,1,2 from left to right.  We start in lane 1
  int target_lane = 1;   //current lane, or lane being moved into.  target_lane != old_lane denotes a lane changing state
  vector<double> d_path(1000,6.0); //d values to be used during a lane change.  Can set to large initial size; will be reset when lane change starts
  h.onMessage([&old_lane,&target_lane,&ref_vel,&d_path,&map_waypoints_x,&map_waypoints_y,&map_waypoints_s,&map_waypoints_dx,&map_waypoints_dy](uWS::WebSocket<uWS::SERVER> ws, char *data, size_t length,uWS::OpCode opCode) {
    // "42" at the start of the message means there's a websocket message event.
    // The 4 signifies a websocket message
    // The 2 signifies a websocket event
    //auto sdata = string(data).substr(0, length);
    //cout << sdata << endl;
    //cout<<"received request from simulator\n";

    if (length && length > 2 && data[0] == '4' && data[1] == '2') {

      auto s = hasData(data);

      if (s != "") {
        auto j = json::parse(s);
        
        string event = j[0].get<string>();
        
        if (event == "telemetry") {
          // j[1] is the data JSON object
          
        	// Main car's localization Data
          	double car_x = j[1]["x"];
          	double car_y = j[1]["y"];
          	double car_s = j[1]["s"];
          	double car_d = j[1]["d"];
          	double car_yaw = j[1]["yaw"];
          	double car_speed = j[1]["speed"];


          	// Previous path data given to the Planner
          	auto previous_path_x = j[1]["previous_path_x"];
          	auto previous_path_y = j[1]["previous_path_y"];
          	// Previous path's end s and d values 
          	double end_path_s = j[1]["end_path_s"];
          	double end_path_d = j[1]["end_path_d"];

          	// Sensor Fusion Data, a list of all other cars on the same side of the road.
          	auto sensor_fusion = j[1]["sensor_fusion"];

          	json msgJson;

                //useful values for computation
                int num_points = 25;             //number of points to send in each response.  Since keeping all previous points, don't want this to be too large and slow reaction time
                double lane_change_time = 3.0;   //time to spend changing lanes in seconds (but note that the time outside of a single lane is smaller than this)
                double lane_width = 4.0;         //lanes are 4 meters wide
                double time_step = 0.02;         //timestep is 0.02 seconds
                double max_vel = 49.5 * 0.44704; // target 49.5 mph, convert to m/s
                double max_accl_step = 0.15;     //max acceleration I'll apply, in m/(s^2)
                int num_spline_points = 1000;    //number of points to generate d values for.  Must be long enough to cover full spline
                //non-persistent values for computation
                int change_dir = target_lane - old_lane;     //direction of lane change (0 if not changing lanes, -1 for left, 1 for right)
                int prev_size = previous_path_x.size();      //number of points remaining from previously sent path
                int steps_advanced = num_points - prev_size; //number of points used by simulator in previously sent path
          	vector<double> next_x_vals;                  //x coordinates of points to be sent
          	vector<double> next_y_vals;                  //y coordinates of points to be sent
                int num_lane_change_points = (int) lane_change_time / time_step; //number of d points to generate from jmt for a lane change.  smaller than spline points

                //cout <<"car_x,car_y are "<<car_x<<","<<car_y<<", ref_vel is "<<ref_vel<<", change_dir is "<<change_dir<<"\n";
                
          	// TODO: define a path made up of (x,y) points that the car will visit sequentially every .02 seconds

                //delete points in d_path that were covered
                d_path.erase(d_path.begin(),d_path.begin()+steps_advanced);
                //fill up d_path with target_lane
                while (d_path.size() < num_spline_points){
                  d_path.push_back(lane_width * (target_lane+0.5));
                }

                //PREDICT LOCATION OF ALL OTHER CARS AT OUR REFERENCE TIME POINT
                double ref_s = car_s;
                if(prev_size > 0){
                  ref_s = end_path_s; //apply predition at refrence point
                }
                double pred_s = ref_s + ref_vel * lane_change_time;      //prediction of where our car will be at the end of lane change time
                double pred2_s = ref_s + ref_vel * 2 * lane_change_time; //prediction of where our car will be at the end of 2xlane change time (for double lane change)
                vector<double>others_id;
                vector<double>others_s;
                vector<double>others_ref_s;
                vector<double>others_pred_s;
                vector<double>others_pred2_s;
                vector<double> others_ref_d;
                vector<double> others_speed;
                for (int i = 0; i < sensor_fusion.size(); i++){ //iterate through all the other cars reported
                  //each sensor_fusion entry is an array of id,x,y,vx,vy,s,d
                  int    other_id = sensor_fusion[i][0];
                  double other_vx = sensor_fusion[i][3];
                  double other_vy = sensor_fusion[i][4];
                  double other_s  = sensor_fusion[i][5];
                  double other_d  = sensor_fusion[i][6];
                  double other_speed = distance(0,0,other_vx,other_vy);
                  others_id.push_back(other_id);
                  others_s.push_back(other_s);
                  double other_ref_s = other_s + ((double)prev_size * time_step * other_speed);//predict where other car will be at our ref point.  Assume car continues in its current lane
                  others_ref_s.push_back(other_ref_s);
                  double other_pred_s = other_ref_s +(lane_change_time * other_speed); //prediction of where other car will be at end of a lane change
                  others_pred_s.push_back(other_pred_s);
                  double other_pred2_s = other_ref_s +(2 * lane_change_time * other_speed); //prediction of where other car will be after 2 lane changes, for double lane change
                  others_pred2_s.push_back(other_pred2_s);
                  others_ref_d.push_back(other_d);  //assume other car isn't moving laterally; this holds reasonably well
                  others_speed.push_back(other_speed);
                }

                //CHECK LANES FOR CLOSE CARS; SET SPEED FOR EACH LANE; MARK IF SAFE TO ENTER LANE
                vector<bool> lane_clear(3,true);
                vector<bool> lane_clear2(3,true);
                vector<double> lane_speed(3,max_vel);
                vector<double> lane_cur_speed(3,max_vel);
                vector<bool> way_too_close(3,false);
                for (int i = 0; i < others_s.size(); i++){
                  double other_d = others_ref_d[i];
                  double other_s = others_s[i];
                  double other_ref_s = others_ref_s[i];
                  double other_pred_s = others_pred_s[i];
                  double other_pred2_s = others_pred2_s[i];
                  double other_speed = others_speed[i];
                  int other_lane = (int)(other_d/lane_width);
                  double other_dist = other_s - car_s;
                  double other_ref_dist = other_ref_s - ref_s;
                  double other_pred_dist = other_pred_s - pred_s;
                  double other_pred2_dist = other_pred2_s - pred2_s;
                  
                  if(other_ref_dist > 0.0 && other_ref_dist < 80.0 && lane_speed[other_lane] > other_speed){
                    //set speed of lane to the slowest car in front of us, within 80m
                    lane_speed[other_lane] = other_speed;
                  }
                  if(other_ref_dist > 0.0 && other_ref_dist < 30.0 && lane_cur_speed[other_lane] > other_speed){
                    //set speed we need to be going in this lane to slowest car in front of us, within 30m
                    lane_cur_speed[other_lane] = other_speed;
                  }
                  if(other_ref_dist > 0.0 && other_ref_dist < 15.0){
                    //flag to slow down hard, if other conditions are met
                    way_too_close[other_lane] = true;
                  }
                  if((other_ref_dist > -15.0 && other_ref_dist < 15.0) || (other_pred_dist > -15.0 && other_pred_dist < 15.0)){
                    //if there is a car within 15m at the reference point, or if there will be at the end of a lane change, we can't move into this lane
                    lane_clear[other_lane] = false;
                  }
                  if((other_pred_dist > -15.0 && other_pred_dist < 15.0) || (other_pred2_dist > -15.0 && other_pred2_dist < 15.0)){
                    //if there will be a car within 15m after 1 lane change, or will be after 2 lane changes, we can't count on a double lane change to this lane
                    lane_clear2[other_lane] = false;
                  }
                }//for i over sensor_fusion

                //SET TARGET VELOCITY BASED ON LANE SPEEDS
                double target_vel = lane_cur_speed[old_lane]; //set target speed to the speed of our current lane
                if (lane_cur_speed[target_lane] < target_vel){ //oops: target lane is now slower than old lane, but I'll finish the maneuver.  Need to slow down
                  target_vel = lane_cur_speed[target_lane];
                  cout <<"    NOTE: target_vel reduced to "<<target_vel<<" due to target_lane "<<target_lane<<" having a slower speed (this should be rare)\n";
                } 

                //DECIDE WHETHER TO CHANGE LANES
                bool initiate_lane_change = false;
                for (int i = 0; i < lane_clear.size(); i++){
                  if(
                     change_dir == 0                      &&           //not already changing lanes
                     lane_clear[i]                        &&           //prospective lane is clear
                     lane_speed[i] > lane_speed[target_lane] + (i == 1 ? -1e-3 : 1.0) &&  //prosepective lane is faster than current lane, and by enough to justify changing- but have a preference for lane 1, since it can access all 3 lanes
                     abs(old_lane - i) < 1.1                           //prospective lane next to current lane
                     ){
                    target_lane = i;
                    initiate_lane_change = true;
                    //cout <<"    lane is "<<old_lane<<"; changed target_lane to "<<target_lane<<"\n";
                  }//meet conditions to change lanes
                }//for i over lane count

                //consider a double lane change: eg, if in lane 0 and lane 1 is slower but lane 2 is faster, and all are clear, move over
                int lane2 = 2 - old_lane;
                if (change_dir == 0 &&  //check to see if 2 lanes over is faster and both lanes are clear; if so, change to lane 1
                    (old_lane == 0 || old_lane == 2) &&
                    lane_clear[1] &&
                    lane_clear2[lane2] &&
                    lane_speed[lane2] > lane_speed[target_lane] + 2.0
                    ){
                  target_lane = 1;
                  initiate_lane_change = true;
                  //cout <<"    lane is "<<old_lane<<"; changed target_lane to "<<target_lane<<" in preparation of a double lane change\n";
                }//if <check for starting double lane change>

                //IF CHANGING LANES, RECOMPUTE d_path 
                if (initiate_lane_change){
                  change_dir = target_lane - old_lane;
                  d_path.erase(d_path.begin()+prev_size,d_path.end()); //erase d_path from before reference point
                  vector<double> d_lane_change = compute_lane_change(old_lane,target_lane,lane_width,num_lane_change_points); //generate d values for jmt
                  d_path.insert(d_path.end(),d_lane_change.begin(),d_lane_change.end()); //concatenate d values for lane change to end of d_path
                }
                    
                //fill up d_lange_change again, in case it was reinitialized, so that it can always be used for spline
                while (d_path.size() < num_spline_points){
                  d_path.push_back(lane_width * (target_lane+0.5));
                }
                
                //CHECK IF LANE CHANGE IS COMPLETE
                if(change_dir != 0){
                  double lane_error = car_d / lane_width - 0.5 - target_lane;
                  if(lane_error > -0.1 && lane_error < 0.1){
                    //cout <<"    lane change completed \n";
                    old_lane = target_lane;
                    change_dir = 0;
                    cout <<"    lane change complete; now in lane "<<old_lane<<"\n";
                  }
                }
                
                //INITIALIZE PATH, OR COPY PREVIOUS PATH
                vector<double> ptsx;
                vector<double> ptsy;
                double ref_x = car_x;
                double ref_y = car_y;
                double ref_x_prev = car_x - 1e-2 * cos(deg2rad(car_yaw));//derive fake previous point from car's heading
                double ref_y_prev = car_y - 1e-2 * sin(deg2rad(car_yaw));
                double ref_yaw = deg2rad(car_yaw);
                if(prev_size >= 2){
                  //use tail end of previous path
                  ref_x = previous_path_x[prev_size-1];
                  ref_y = previous_path_y[prev_size-1];
                  double ref_x_prev = previous_path_x[prev_size - 2];
                  double ref_y_prev = previous_path_y[prev_size - 2];
                  ref_yaw = atan2(ref_y-ref_y_prev,ref_x-ref_x_prev);
                }
                ptsx.push_back(ref_x_prev);
                ptsx.push_back(ref_x);
                ptsy.push_back(ref_y_prev);
                ptsy.push_back(ref_y);

                //COMPUTE NEW PATH
                for (int i = 1; i <= 3; i++){
                  //first step must exceed end of previous path; compare dist_ahead to max_vel and num_points
                  double dist_ahead = i * 40.0;
                  int num_points_ahead;
                  if (ref_vel > 10.0){
                    num_points_ahead = (int) (dist_ahead / ref_vel / time_step);
                  }else{
                    num_points_ahead = i * 100; //won't matter much since we're going so slow
                  }
                  double d_ahead = d_path[num_points_ahead];
                  //cout <<"    creating spline point s,d = "<<dist_ahead<<","<<d_ahead<<"\n";
                  vector<double> next_wp = getXY(car_s+dist_ahead,d_ahead,map_waypoints_s,map_waypoints_x,map_waypoints_y);
                  ptsx.push_back(next_wp[0]);
                  ptsy.push_back(next_wp[1]);

                }
                //convert points to xy points in the car's frame of reference
                for (int i=0; i<ptsx.size();i++){
                  double shift_x = ptsx[i]-ref_x;
                  double shift_y = ptsy[i]-ref_y;
                  ptsx[i] = shift_x*cos(0.0-ref_yaw) - shift_y*sin(0.0-ref_yaw);
                  ptsy[i] = shift_x*sin(0.0-ref_yaw) + shift_y*cos(0.0-ref_yaw);
                }

                tk::spline spline;
                spline.set_points(ptsx,ptsy);

                //copy previous points to new path
                for(int i = 0; i < prev_size; i++){
                  next_x_vals.push_back(previous_path_x[i]);
                  next_y_vals.push_back(previous_path_y[i]);
                }
                double target_x = 30.0; //meters out to look when adjusting distance for spline
                double target_y = spline(target_x);
                double target_dist = distance(0.0,0.0,target_x,target_y);
                double target_fact = target_dist / target_x; //factor to adjust distance in x direction to aproximate spline distance in 2d space
                double ref_accl_step = 0.0;                  //acceleration to apply at reference point

                //compute what our acceleration should be, based on way_too_close, ref_vel, and target_vel
                if (way_too_close[old_lane] || way_too_close[target_lane]) { 
                  ref_accl_step = -1 * max_accl_step;
                }else{
                  ref_accl_step = max_accl_step * (target_vel - ref_vel) / (1.0 * target_vel); //scale acceleration based on difference in speed gradually
                  if (fabs(ref_accl_step) > max_accl_step){
                    ref_accl_step = max_accl_step * (ref_accl_step > 0 ? 1 : -1);
                  }
                }

                double iter_x = 0.0;
                for (int i = 1; i <= num_points - prev_size; i++){

                  ref_vel += ref_accl_step; //increase velocity at every step, for constant acceleration
                  double ref_vel_dist = time_step * ref_vel;
                  //cout<<"    set ref_vel,ref_vel_dist to "<<ref_vel<<","<<ref_vel_dist<<" since target_vel,way_too_close are "<<target_vel<<","<<way_too_close[old_lane]<<","<<way_too_close[target_lane]<<"\n";

                  iter_x += ref_vel_dist/target_fact; //coordinate conversion set ref_x to 0, which is where we want to add points to
                  double iter_y = spline(iter_x);

                  //convert point back to global xy coordinates
                  double x_point = iter_x*cos(ref_yaw) - iter_y*sin(ref_yaw) + ref_x;
                  double y_point = iter_x*sin(ref_yaw) + iter_y*cos(ref_yaw) + ref_y;

                  next_x_vals.push_back(x_point);
                  next_y_vals.push_back(y_point);
                }

                //set the message to our computed points!
                msgJson["next_x"] = next_x_vals;
          	msgJson["next_y"] = next_y_vals;

          	auto msg = "42[\"control\","+ msgJson.dump()+"]";

          	//this_thread::sleep_for(chrono::milliseconds(1000));
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
  if (h.listen(port)) {
    std::cout << "Listening to port " << port << std::endl;
  } else {
    std::cerr << "Failed to listen to port" << std::endl;
    return -1;
  }
  h.run();
}
