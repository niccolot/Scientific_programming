#include <stdio.h>
#include <math.h>
#include <stdbool.h>

/*
script for studying the behaviour of different ODE integration algorithms 
for 1d second order differential equations using the harmonic oscillator
as example x''(t) = -omega2*x(t)

the ode can be modified just by changing 'get_force' and 'get_energy' functions

selecting the boolean variables 'TRAJECTORIES' and 'INT_STEP_DEPENDECE' one
can choose whether to study the trajectories in function of t at fixed integration step dt
or the dependence of the reduce energy (E(T)-E(0))/E(0) in function of dt

writes in files the trajectories (t,x(t),v(t)), the energies (t,E(t)) and reduce energies
(dt, (E(T)-E(0))/E(0))
*/

struct PhaseSpace{
    double x;
    double v;
};


struct PhaseSpace *Euler_init(double x0, double v0);
struct PhaseSpace *Euler_Cromer_init(double x0, double v0);
struct PhaseSpace *Velocity_Verlet_init(double x0, double v0);
struct PhaseSpace *Runge_Kutta_init(double x0, double v0);
struct PhaseSpace *Euler(double dt, double omega2, struct PhaseSpace *phase_space_old);
struct PhaseSpace *Euler_Cromer(double dt, double omega2, struct PhaseSpace *phase_space_old);
struct PhaseSpace *Velocity_Verlet(double dt, double omega2, struct PhaseSpace *phase_space_old);
struct PhaseSpace *Runge_Kutta(double dt, double omega2, struct PhaseSpace *phase_space_old); // second order rk
double get_energy(double omega2, struct PhaseSpace *phase_space);
double get_force(double omega2, double x);


int main(){

    bool TRAJECTORIES = true;
    bool INT_STEP_DEPENDENCE = true;
    double dt = 0.01;
    double dt_list[] = {0.001, 0.002, 0.003, 0.005, 0.01, 0.015, 0.02, 0.05, 0.075, 0.1}; 
    int T0 = 0;
    long int T = 10000; // integration steps, not seconds
    long int T_list[] = {1000, 500, 333, 200, 100, 66, 50, 20, 13, 10}; // 1 second (T*dt) of simulation for each dt
    double omega2 = 0.5;
    double x0 = 1.0;
    double v0 = 0.0;
    double e0;
    double energy_euler;
    double energy_euler_cromer;
    double energy_velocity_verlet;
    double energy_runge_kutta;
    struct PhaseSpace euler;
    struct PhaseSpace euler_cromer;
    struct PhaseSpace velocity_verlet;
    struct PhaseSpace runge_kutta;
    // trajectories
    FILE *pf_trajectory_euler;
    FILE *pf_trajectory_euler_cromer;
    FILE *pf_trajectory_vv;
    FILE *pf_trajectory_rk;
    // energy in function of t, E(t)
    FILE *pf_energy_vv;
    FILE *pf_energy_rk;
    FILE *pf_energy_euler;
    FILE *pf_energy_euler_cromer;
    // energy at fixed T varying the integration step
    FILE *pf_energy_dt_e;
    FILE *pf_energy_dt_ec;
    FILE *pf_energy_dt_vv;
    FILE *pf_energy_dt_rk;


    if(TRAJECTORIES){

        pf_trajectory_euler = fopen("trajectory_euler.txt", "w");
        pf_trajectory_euler_cromer = fopen("trajectory_euler_cromer.txt", "w");
        pf_trajectory_vv = fopen("trajectory_vv.txt", "w");
        pf_trajectory_rk = fopen("trajectory_rk.txt", "w");
        pf_energy_euler = fopen("energy_euler.txt", "w");
        pf_energy_euler_cromer = fopen("energy_euler_cromer.txt", "w");
        pf_energy_vv = fopen("energy_vv.txt", "w");
        pf_energy_rk = fopen("energy_rk.txt", "w");

        euler = *Euler_Cromer_init(x0,v0);
        euler_cromer = *Euler_Cromer_init(x0,v0);
        velocity_verlet = *Velocity_Verlet_init(x0,v0);
        runge_kutta = *Runge_Kutta_init(x0,v0);

        for(int t=T0;t<T;t++){

            euler = *Euler(dt, omega2, &euler);
            euler_cromer = *Euler_Cromer(dt, omega2, &euler_cromer);
            velocity_verlet = *Velocity_Verlet(dt, omega2, &velocity_verlet);
            runge_kutta = *Runge_Kutta(dt, omega2, &runge_kutta);

            energy_euler = get_energy(omega2, &euler);
            energy_euler_cromer = get_energy(omega2, &euler_cromer);
            energy_velocity_verlet = get_energy(omega2, &velocity_verlet);
            energy_runge_kutta = get_energy(omega2, &runge_kutta);

            fprintf(pf_trajectory_euler, "%d\t%f\t%f\n", t, euler.x, euler.v);
            fprintf(pf_trajectory_euler_cromer, "%d\t%f\t%f\n", t, euler_cromer.x, euler_cromer.v);
            fprintf(pf_trajectory_vv, "%d\t%f\t%f\n", t, velocity_verlet.x, velocity_verlet.v);
            fprintf(pf_trajectory_rk, "%d\t%f\t%f\n", t, runge_kutta.x, runge_kutta.v);
            fprintf(pf_energy_euler, "%d\t%f\n", t, energy_euler);
            fprintf(pf_energy_euler_cromer, "%d\t%f\n", t, energy_euler_cromer);
            fprintf(pf_energy_vv, "%d\t%f\n", t, energy_velocity_verlet);
            fprintf(pf_energy_rk, "%d\t%f\n", t, energy_runge_kutta);
        }

        fclose(pf_trajectory_euler);
        fclose(pf_trajectory_euler_cromer);
        fclose(pf_trajectory_vv);
        fclose(pf_trajectory_rk);
        fclose(pf_energy_euler);
        fclose(pf_energy_euler_cromer);
        fclose(pf_energy_vv);
        fclose(pf_energy_rk);
    }
    
    if(INT_STEP_DEPENDENCE){

        pf_energy_dt_e = fopen("energy_dt_e.txt", "w");
        pf_energy_dt_ec = fopen("energy_dt_ec.txt", "w");
        pf_energy_dt_vv = fopen("energy_dt_vv.txt", "w");
        pf_energy_dt_rk = fopen("energy_dt_rk.txt", "w");

        euler = *Euler_Cromer_init(x0,v0);
        euler_cromer = *Euler_Cromer_init(x0,v0);
        velocity_verlet = *Velocity_Verlet_init(x0,v0);
        runge_kutta = *Runge_Kutta_init(x0,v0);

        e0 = get_energy(omega2, &euler); //initial energy is the same for all algorithms

        // 10 different integration steps 'dt' used in dt_lsit
        for(int i=0;i<10;i++){

            dt = dt_list[i];
            T = T_list[i];
            
            for(int t=0;t<T;t++){

                euler = *Euler(dt, omega2, &euler);
                euler_cromer = *Euler_Cromer(dt, omega2, &euler_cromer);
                velocity_verlet = *Velocity_Verlet(dt, omega2, &velocity_verlet);
                runge_kutta = *Runge_Kutta(dt, omega2, &runge_kutta);
            }
            
            energy_euler = get_energy(omega2, &euler);
            energy_euler_cromer = get_energy(omega2, &euler_cromer);
            energy_velocity_verlet = get_energy(omega2, &velocity_verlet);
            energy_runge_kutta = get_energy(omega2, &runge_kutta);

            energy_euler = (energy_euler-e0)/e0;
            energy_euler_cromer = (energy_euler_cromer-e0)/e0;;
            energy_velocity_verlet = (energy_velocity_verlet-e0)/e0;;
            energy_runge_kutta = (energy_runge_kutta-e0)/e0;;

            fprintf(pf_energy_dt_e, "%f\t%f\n", dt, energy_euler);
            fprintf(pf_energy_dt_ec, "%f\t%f\n", dt, energy_euler_cromer);
            fprintf(pf_energy_dt_vv, "%f\t%f\n", dt, energy_velocity_verlet);
            fprintf(pf_energy_dt_rk, "%f\t%f\n", dt, energy_runge_kutta);
        }

        fclose(pf_energy_dt_e);
        fclose(pf_energy_dt_ec);
        fclose(pf_energy_dt_vv);
        fclose(pf_energy_dt_rk);
    }   
    
    return 0;
}


////////////////////////////////////FUNCTIONS////////////////////////////////////////////////
struct PhaseSpace *Euler_init(double x0, double v0){

    static struct PhaseSpace euler_init;

    euler_init.x = x0;
    euler_init.v = v0;

    return &euler_init;
}


struct PhaseSpace *Euler_Cromer_init(double x0, double v0){

    static struct PhaseSpace euler_cromer_init;

    euler_cromer_init.x = x0;
    euler_cromer_init.v = v0;

    return &euler_cromer_init;
}


struct PhaseSpace *Velocity_Verlet_init(double x0, double v0){

    static struct PhaseSpace velocity_verlet_init;

    velocity_verlet_init.x = x0;
    velocity_verlet_init.v = v0;

    return &velocity_verlet_init;
}


struct PhaseSpace *Runge_Kutta_init(double x0, double v0){

    static struct PhaseSpace runge_kutta_init;

    runge_kutta_init.x = x0;
    runge_kutta_init.v = v0;

    return &runge_kutta_init;
}


struct PhaseSpace *Euler(double dt, double omega2, struct PhaseSpace *phase_space_old){

    static struct PhaseSpace euler_step;

    euler_step.v = (phase_space_old->v) + dt*get_force(omega2, phase_space_old->x);
    euler_step.x = (phase_space_old->x) + (phase_space_old->v)*dt;

    return &euler_step;
}


struct PhaseSpace *Euler_Cromer(double dt, double omega2, struct PhaseSpace *phase_space_old){

    static struct PhaseSpace euler_cromer_step;

    euler_cromer_step.v = (phase_space_old->v) + dt*get_force(omega2, phase_space_old->x);
    euler_cromer_step.x = (phase_space_old->x) + euler_cromer_step.v*dt;

    return &euler_cromer_step;
}


struct PhaseSpace *Velocity_Verlet(double dt, double omega2, struct PhaseSpace *phase_space_old){

    static struct PhaseSpace velocity_verlet_step;
    double f_old = get_force(omega2, phase_space_old->x);

    velocity_verlet_step.x = (phase_space_old->x) + (phase_space_old->v)*dt + 0.5*f_old*pow(dt,2);
    
    double f_new = get_force(omega2, velocity_verlet_step.x);

    velocity_verlet_step.v = (phase_space_old->v) + 0.5*(f_old+f_new)*dt;

    return &velocity_verlet_step;
}


struct PhaseSpace *Runge_Kutta(double dt, double omega2, struct PhaseSpace *phase_space_old){

    static struct PhaseSpace runge_kutta_step;
    double x_p = (phase_space_old->v)*dt;
    double v_p = get_force(omega2, phase_space_old->x)*dt;
    double x_pp = phase_space_old->x + 0.5*x_p;
    double v_pp = phase_space_old->v + 0.5*v_p;

    runge_kutta_step.x = (phase_space_old->x) + v_pp*dt;
    runge_kutta_step.v = (phase_space_old->v) + get_force(omega2, x_pp)*dt;

    return &runge_kutta_step;
}


double get_energy(double omega2, struct PhaseSpace *phase_space){
    /*
    E(t) = 1/2*(v(t)**2) + 1/2*omega2*(x(t)**2) for the harrmonic oscillator
    */
    double energy = 0.5*pow((phase_space->v),2) + 0.5*omega2*pow((phase_space->x),2);

    return energy;
}


double get_force(double omega2, double x){
    /*
    returns harmonic force for the harmonic oscillator
    */
    double f = -omega2*x;

    return f;
}