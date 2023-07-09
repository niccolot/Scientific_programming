/*
numerical integration for a system corresponding to 2 bodies interacting with a coulomb type interaction
both immersed in a central force filed, e.g. two planets going around with a fixed sun
*/

#include <stdio.h>
#include <math.h>

// struct for the 2d phase space corresponding to 1 degree of freedom of one body
struct PhaseSpace{
    double x;
    double v;
};

// struct for the 4d phase space corresponding to both degree of freedom of 1 body
struct OneBody{
    struct PhaseSpace coo1;
    struct PhaseSpace coo2;
};

// struct for the total 8d phase space
struct TwoBodies{
    struct OneBody body1;
    struct OneBody body2;
};


double get_force1(double mb, double x1, double y1, double x2, double y2);
double get_force2(double mb, double x1, double y1, double x2, double y2);
double get_force3(double ma, double x1, double y1, double x2, double y2);
double get_force4(double ma, double x1, double y1, double x2, double y2);
struct TwoBodies *Runge_Kutta_init(double body1_init[], double body2_init[]);
struct TwoBodies *Runge_Kutta_step(double dt, double ma, double mb, struct TwoBodies *old_system);

int main(){

    // initial conditions -> {x0,vx0,y0,vy0}
    double body1_init[] = {1,0,2,0.1};
    double body2_init[] = {2,0,2,0};
    double dt = 0.001;
    double T = 10000;
    // masses
    double ma = 0.001;
    double mb = 0.01;
    struct TwoBodies system;
    FILE *pf_trajectories1; // body 1
    FILE *pf_trajectories2; // body 2

    pf_trajectories1 = fopen("trajectories1.txt", "w");
    pf_trajectories2 = fopen("trajectories2.txt", "w");

    system = *Runge_Kutta_init(body1_init, body2_init);

    double x1;
    double vx1;
    double y1;
    double vy1;
    double x2;
    double vx2;
    double y2;
    double vy2;

    for(int t=0;t<T;t++){

        // integration step
        system = *Runge_Kutta_step(dt, ma, mb, &system);
        
        // new states
        x1 = system.body1.coo1.x;
        vx1 = system.body1.coo1.v;
        y1 = system.body1.coo2.x;
        vy1 = system.body1.coo2.v;
        x2 = system.body2.coo1.x;
        vx2 = system.body2.coo1.v;
        y2 = system.body2.coo2.x;
        vy2 = system.body2.coo2.v;

        // save the new states as (t,x(t),vx(t),y(t),vy(t))
        fprintf(pf_trajectories1, "%d\t%f\t%f\t%f\t%f\n", t,x1,vx1,y1,vy1);
        fprintf(pf_trajectories2, "%d\t%f\t%f\t%f\t%f\n", t,x2,vx2,y2,vy2);
    }

    fclose(pf_trajectories1);
    fclose(pf_trajectories2);

    return 0;
}


struct TwoBodies *Runge_Kutta_init(double body1_init[], double body2_init[]){

    static struct TwoBodies init;

    init.body1.coo1.x = body1_init[0];
    init.body1.coo1.v = body1_init[1];
    init.body1.coo2.x = body1_init[2];
    init.body1.coo2.v = body1_init[3];
    init.body2.coo1.x = body2_init[0];
    init.body2.coo1.v = body2_init[1];
    init.body2.coo2.x = body2_init[2];
    init.body2.coo2.v = body2_init[3];

    return &init;
}


struct TwoBodies *Runge_Kutta_step(double dt, double ma, double mb, struct TwoBodies *old_system){
    /*
    second order runge kutta integration step
    */
    static struct TwoBodies rk_step;

    // positions of body 1 and 2
    double x1 = old_system->body1.coo1.x;
    double y1 = old_system->body1.coo2.x;
    double x2 = old_system->body2.coo1.x;
    double y2 = old_system->body2.coo2.x;

    // velocities of body 1 and 2
    double vx1 = old_system->body1.coo1.v;
    double vy1 = old_system->body1.coo2.v;
    double vx2 = old_system->body2.coo1.v;
    double vy2 = old_system->body2.coo2.v;

    // two 2d newtons equations -> 4 forces
    double f1 = get_force1(mb, x1, y1, x2, y2);
    double f2 = get_force1(mb, x1, y1, x2, y2);
    double f3 = get_force1(ma, x1, y1, x2, y2);
    double f4 = get_force1(ma, x1, y1, x2, y2);
    
    // runge kutta parameters ofr body 1 coordinate x
    double x_p1 = vx1*dt;
    double vx_p1 = f1*dt;
    double x_pp1 = x1 + 0.5*x_p1;
    double vx_pp1 = vx1 + 0.5*vx_p1;

    // runge kutta parameters for body 1 coordinate y
    double y_p1 = vy1*dt;
    double vy_p1 = f2*dt;
    double y_pp1 = y1 + 0.5*y_p1;
    double vy_pp1 = vy1 + 0.5*vy_p1;

    // runge kutta parameters for body 2 coordinate x
    double x_p2 = vx2*dt;
    double vx_p2 = f3*dt;
    double x_pp2 = x2 + 0.5*x_p2;
    double vx_pp2 = vx2 + 0.5*vx_p2;

    // runge kutta parameters for body 2 coordinate y
    double y_p2 = vy2*dt;
    double vy_p2 = f4*dt;
    double y_pp2 = y2 + 0.5*y_p2;
    double vy_pp2 = vy2 + 0.5*vy_p2;
    

    // runge kutta algorithm for body 1 coordinate x
    rk_step.body1.coo1.x = x1 + vx_pp1*dt;
    rk_step.body1.coo1.v = vx1 + get_force1(mb,x_pp1,y1,x2,y2)*dt;

    // runge kutta algorithm for body 1 coordinate y
    rk_step.body1.coo2.x = y1 + vy_pp1*dt;
    rk_step.body1.coo2.v = vy1 + get_force2(mb,x1,y_pp1,x2,y2)*dt;

    // runge kutta algorithm for body 2 coordinate x
    rk_step.body2.coo1.x = x2 + vx_pp2*dt;
    rk_step.body2.coo1.v = vx2 + get_force3(ma,x1,y1,x_pp2,y2)*dt;

    // runge kutta algorithm for body 2 coordinate y
    rk_step.body2.coo2.x = y2 + vy_pp2*dt;
    rk_step.body2.coo2.v = vy2 + get_force4(ma,x1,y1,x2,y_pp2)*dt;


    return &rk_step;
}


double get_force1(double mb, double x1, double y1, double x2, double y2){

    double f;
    double r_a = pow(pow(x1,2)+pow(y1,2), 1.5);
    double r_ab = pow(pow(x2-x1, 2)+pow(y2-y1, 2), 1.5);

    f = -x1/r_a + (mb*(x2-x1))/r_ab;
    
    return f;
}


double get_force2(double mb, double x1, double y1, double x2, double y2){

    double f;
    double r_a = pow(pow(x1,2)+pow(y1,2), 1.5);
    double r_ab = pow(pow(x2-x1, 2)+pow(y2-y1, 2), 1.5);

    f = -y1/r_a + (mb*(y2-y1))/r_ab;
    
    return f;
}


double get_force3(double ma, double x1, double y1, double x2, double y2){

    double f;
    double r_b = pow(pow(x2,2)+pow(y2,2), 1.5);
    double r_ab = pow(pow(x2-x1, 2)+pow(y2-y1, 2), 1.5);

    f = -x2/r_b - (ma*(x2-x1))/r_ab;
    
    return f;
}


double get_force4(double ma, double x1, double y1, double x2, double y2){

    double f;
    double r_b = pow(pow(x2,2)+pow(y2,2), 1.5);
    double r_ab = pow(pow(x2-x1, 2)+pow(y2-y1, 2), 1.5);

    f = -y1/r_b - (ma*(y2-y1))/r_ab;
    
    return f;
}