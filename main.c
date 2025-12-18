#include <stdio.h>
#include <math.h>
#define RAYGUI_IMPLEMENTATION
#include "raygui.h"

int mainmenu = 1;
int integrator;
int play = 1;
const double G = 0.2;
const double energyFactor = 6.6e33;

typedef struct
{
    double x, y, vx, vy, m, r;
} Body;

typedef struct
{
    double x, y;
} Force;

Force computeForce(Body a, Body b)
{
    double dx = b.x - a.x;
    double dy = b.y - a.y;
    double r = sqrt(dx * dx + dy * dy + 1e-6);
    double mag = (G * a.m * b.m) / (r * r * r);
    Force f;
    f.x = mag * dx;
    f.y = mag * dy;
    return f;
}

void applyRK2(double dt, Body p[], Force f[])
{
    int n = 9;
    Body temp[9];
    Force f_mid[9];

    // Step 1: Compute mid-step positions and velocities
    for (int i = 0; i < n; i++)
    {
        double ax = f[i].x / p[i].m;
        double ay = f[i].y / p[i].m;

        temp[i].x = p[i].x + p[i].vx * dt * 0.5;
        temp[i].y = p[i].y + p[i].vy * dt * 0.5;
        temp[i].vx = p[i].vx + ax * dt * 0.5;
        temp[i].vy = p[i].vy + ay * dt * 0.5;
        temp[i].m = p[i].m;
    }

    // Step 2: Compute forces at mid-step
    for (int i = 0; i < n; i++)
        f_mid[i] = (Force){0, 0};
    for (int i = 0; i < n; i++)
        for (int j = i + 1; j < n; j++)
        {
            Force fij = computeForce(temp[i], temp[j]);
            f_mid[i].x += fij.x;
            f_mid[i].y += fij.y;
            f_mid[j].x -= fij.x;
            f_mid[j].y -= fij.y;
        }

    // Step 3: Update positions and velocities using mid-step accelerations
    for (int i = 0; i < n; i++)
    {
        double ax_mid = f_mid[i].x / p[i].m;
        double ay_mid = f_mid[i].y / p[i].m;

        // Positions updated using **mid-step velocities**
        p[i].x += (p[i].vx + 0.5 * ax_mid * dt) * dt;
        p[i].y += (p[i].vy + 0.5 * ay_mid * dt) * dt;

        // Velocities updated using **mid-step accelerations**
        p[i].vx += ax_mid * dt;
        p[i].vy += ay_mid * dt;
    }
}

void updatePhysics(int i, int integrator, Force f, double dt, Body p[], Force fo[])
{
    switch (integrator)
    {
    case 1: // Euler
        p[i].x += p[i].vx * dt;
        p[i].y += p[i].vy * dt;
        p[i].vx += (f.x / p[i].m) * dt;
        p[i].vy += (f.y / p[i].m) * dt;
        break;
    case 2: // Symplectic Euler
        p[i].vx += (f.x / p[i].m) * dt;
        p[i].vy += (f.y / p[i].m) * dt;
        p[i].x += p[i].vx * dt;
        p[i].y += p[i].vy * dt;
        break;
    case 4: // RK4
    {
        double x0 = p[i].x, y0 = p[i].y;
        double vx0 = p[i].vx, vy0 = p[i].vy;
        double m = p[i].m;

        double k1_x = vx0, k1_y = vy0;
        double k1_vx = f.x / m, k1_vy = f.y / m;

        Body s2 = {x0 + k1_x * dt / 2, y0 + k1_y * dt / 2, vx0 + k1_vx * dt / 2, vy0 + k1_vy * dt / 2, m, 0};
        Force f2 = {0, 0};
        for (int j = 0; j < 9; j++)
            if (j != i)
            {
                Force ff = computeForce(s2, p[j]);
                f2.x += ff.x;
                f2.y += ff.y;
            }

        double k2_x = s2.vx, k2_y = s2.vy;
        double k2_vx = f2.x / m, k2_vy = f2.y / m;

        Body s3 = {x0 + k2_x * dt / 2, y0 + k2_y * dt / 2, vx0 + k2_vx * dt / 2, vy0 + k2_vy * dt / 2, m, 0};
        Force f3 = {0, 0};
        for (int j = 0; j < 9; j++)
            if (j != i)
            {
                Force ff = computeForce(s3, p[j]);
                f3.x += ff.x;
                f3.y += ff.y;
            }

        double k3_x = s3.vx, k3_y = s3.vy;
        double k3_vx = f3.x / m, k3_vy = f3.y / m;

        Body s4 = {x0 + k3_x * dt, y0 + k3_y * dt, vx0 + k3_vx * dt, vy0 + k3_vy * dt, m, 0};
        Force f4 = {0, 0};
        for (int j = 0; j < 9; j++)
            if (j != i)
            {
                Force ff = computeForce(s4, p[j]);
                f4.x += ff.x;
                f4.y += ff.y;
            }

        double k4_x = s4.vx, k4_y = s4.vy;
        double k4_vx = f4.x / m, k4_vy = f4.y / m;

        p[i].x += dt / 6 * (k1_x + 2 * k2_x + 2 * k3_x + k4_x);
        p[i].y += dt / 6 * (k1_y + 2 * k2_y + 2 * k3_y + k4_y);
        p[i].vx += dt / 6 * (k1_vx + 2 * k2_vx + 2 * k3_vx + k4_vx);
        p[i].vy += dt / 6 * (k1_vy + 2 * k2_vy + 2 * k3_vy + k4_vy);
        break;
    }
    default:
        break;
    }
}

void applyLeapfrog(double dt, Body p[], Force f[])
{
    int n = 9;
    Force newF[9];

    // Half-step velocity (kick)
    for (int i = 0; i < n; i++)
    {
        p[i].vx += f[i].x / p[i].m * dt * 0.5;
        p[i].vy += f[i].y / p[i].m * dt * 0.5;
    }

    // Full-step position (drift)
    for (int i = 0; i < n; i++)
    {
        p[i].x += p[i].vx * dt;
        p[i].y += p[i].vy * dt;
    }

    // Compute forces at new positions
    for (int i = 0; i < n; i++)
        newF[i] = (Force){0, 0};
    for (int i = 0; i < n; i++)
        for (int j = i + 1; j < n; j++)
        {
            Force fij = computeForce(p[i], p[j]);
            newF[i].x += fij.x;
            newF[i].y += fij.y;
            newF[j].x -= fij.x;
            newF[j].y -= fij.y;
        }

    // Half-step velocity (kick)
    for (int i = 0; i < n; i++)
    {
        p[i].vx += newF[i].x / p[i].m * dt * 0.5;
        p[i].vy += newF[i].y / p[i].m * dt * 0.5;
    }

    // Copy new forces back
    for (int i = 0; i < n; i++)
        f[i] = newF[i];
}

void initialize(Body planet[])
{
    int sun_index = 8;
    planet[0].x = 648;
    planet[0].y = 360;
    planet[0].vx = 0;
    planet[0].vy = 0;
    planet[0].m = 0.0000165;
    planet[0].r = 0;
    planet[1].x = 654;
    planet[1].y = 360;
    planet[1].vx = 0;
    planet[1].vy = 0;
    planet[1].m = 0.000245;
    planet[1].r = 0;
    planet[2].x = 660;
    planet[2].y = 360;
    planet[2].vx = 0;
    planet[2].vy = 0;
    planet[2].m = 0.0003;
    planet[2].r = 0;
    planet[3].x = 670;
    planet[3].y = 360;
    planet[3].vx = 0;
    planet[3].vy = 0;
    planet[3].m = 0.000032;
    planet[3].r = 0;
    planet[4].x = 744;
    planet[4].y = 360;
    planet[4].vx = 0;
    planet[4].vy = 0;
    planet[4].m = 0.0954;
    planet[4].r = 0;
    planet[5].x = 832;
    planet[5].y = 360;
    planet[5].vx = 0;
    planet[5].vy = 0;
    planet[5].m = 0.0544;
    planet[5].r = 0;
    planet[6].x = 1024;
    planet[6].y = 360;
    planet[6].vx = 0;
    planet[6].vy = 0;
    planet[6].m = 0.0081;
    planet[6].r = 0;
    planet[7].x = 1241;
    planet[7].y = 360;
    planet[7].vx = 0;
    planet[7].vy = 0;
    planet[7].m = 0.0102;
    planet[7].r = 0;
    planet[8].x = 640;
    planet[8].y = 360;
    planet[8].vx = 0;
    planet[8].vy = 0;
    planet[8].m = 100;
    planet[8].r = 0;

    // Initialize circular velocities around Sun
    for (int i = 0; i < 8; i++)
    {
        double dx = planet[i].x - planet[sun_index].x;
        double r = fabs(dx);
        double v = sqrt(G * planet[sun_index].m / r);
        planet[i].vy = (dx < 0) ? v : -v;
    }

    // Fix Sun momentum so total momentum = 0
    double Px = 0, Py = 0;
    for (int i = 0; i < 8; i++)
    {
        Px += planet[i].m * planet[i].vx;
        Py += planet[i].m * planet[i].vy;
    }
    planet[sun_index].vx = -Px / planet[sun_index].m;
    planet[sun_index].vy = -Py / planet[sun_index].m;
}

double KE(Body p[], int k)
{
    return 0.5 * p[k].m * (p[k].vx * p[k].vx + p[k].vy * p[k].vy) * energyFactor;
}

double PE(Body p[], int k)
{
    double u = 0;
    for (int i = 0; i < 9; i++)
    {
        if (i != k)
        {
            double dx = p[k].x - p[i].x, dy = p[k].y - p[i].y;
            double r = sqrt(dx * dx + dy * dy + 1e-6);
            u -= G * p[k].m * p[i].m * energyFactor / r;
        }
    }
    return u;
}

// ----------------- Main function and GUI remain unchanged -----------------
int main()
{
    InitWindow(1280, 720, "Solar System");
    InitAudioDevice();
    Sound clickSound = LoadSound("/home/manas/Downloads/click-6.mp3");
    Font myFont = LoadFont("/home/manas/Desktop/Body simulator/primus/Primus Bold.otf");
    SetTargetFPS(190);
    Body planet[9];
    initialize(planet);
    Force force[9];
    double dt;
    GuiSetFont(myFont);
    GuiSetStyle(DEFAULT, TEXT_SIZE, 15);
    double TE = 0, TKE = 0, TU = 0;
    FILE *fp = fopen("SIMULATIONdata.csv", "w");
    fprintf(fp, "i,time,x,y,vx,vy,KE,PE,TE\n");
    double t=0;
    while (!WindowShouldClose())
    {
        
        TE = TKE = TU = 0;
        dt = 0.002;
        BeginDrawing();
        ClearBackground(RAYWHITE);

        // Compute energies
        for (int i = 0; i < 9; i++)
            TKE += KE(planet, i);
        for (int i = 0; i < 9; i++)
            for (int j = i + 1; j < 9; j++)
            {
                double dx = planet[j].x - planet[i].x;
                double dy = planet[j].y - planet[i].y;
                double r = sqrt(dx * dx + dy * dy + 1e-6);
                TU -= G * planet[i].m * planet[j].m * energyFactor / r;
            }
        TE = TKE + TU;

        // GUI menu and buttons
        if (mainmenu)
        {
            DrawText("Please Select Numerical Integrator", 500, 220, 19, WHITE);
            if (GuiButton((Rectangle){610, 250, 80, 30}, "Euler"))
            {
                PlaySound(clickSound);
                mainmenu = 0;
                integrator = 1;
            }
            else if (GuiButton((Rectangle){565, 290, 170, 30}, "Simplectic Euler"))
            {
                PlaySound(clickSound);
                mainmenu = 0;
                integrator = 2;
            }
            else if (GuiButton((Rectangle){610, 330, 80, 30}, "RK2"))
            {
                PlaySound(clickSound);
                mainmenu = 0;
                integrator = 3;
            }
            else if (GuiButton((Rectangle){610, 370, 80, 30}, "RK4"))
            {
                PlaySound(clickSound);
                mainmenu = 0;
                integrator = 4;
            }
            else if (GuiButton((Rectangle){610, 410, 80, 30}, "Leepfrog"))
            {
                PlaySound(clickSound);
                mainmenu = 0;
                integrator = 5;
            }
        }
        else
        {
            DrawText(TextFormat("Energy: %.8e", TE), 30, 140, 16, BLACK);
            if (GuiButton((Rectangle){13, 10, 175, 30}, "Change Integrator"))
            {
                mainmenu = 1;
                PlaySound(clickSound);
                initialize(planet);
            }
            if (GuiButton((Rectangle){13, 45, 175, 30}, "Reset Simulation"))
            {
                initialize(planet);
                PlaySound(clickSound);
            }
            if (!play)
            {
                if (GuiButton((Rectangle){13, 80, 175, 30}, "Play Simulation"))
                {
                    play = 1;
                    PlaySound(clickSound);
                }
            }
            else
            {
                if (GuiButton((Rectangle){13, 80, 175, 30}, "Pause Simulation"))
                {
                    play = 0;
                    PlaySound(clickSound);
                }
            }

            if (play)
            {
                for(int i = 0; i < 9; i++)
    {
        fprintf(fp, "%d,%.6f,%.6f,%.6f,%.6f,%.6f,%.6e,%.6e,%.6e\n",
            i, t,
            planet[i].x, planet[i].y,
            planet[i].vx, planet[i].vy,
            KE(planet, i),
            PE(planet, i),
            KE(planet, i) + PE(planet, i)
        );
    }

    t += dt;
                for (int i = 0; i < 9; i++)
                    force[i] = (Force){0, 0};
                for (int i = 0; i < 9; i++)
                    for (int j = i + 1; j < 9; j++)
                    {
                        Force f = computeForce(planet[i], planet[j]);
                        force[i].x += f.x;
                        force[i].y += f.y;
                        force[j].x -= f.x;
                        force[j].y -= f.y;
                    }

                if (integrator == 3)
                    applyRK2(dt, planet, force);
                else if (integrator == 5)
                    applyLeapfrog(dt, planet, force);
                else
                    for (int i = 0; i < 9; i++)
                        updatePhysics(i, integrator, force[i], dt, planet, force);
            }

            for (int i = 0; i < 9; i++)
                DrawCircle(planet[i].x, planet[i].y, 2, GRAY);
        }
        

        EndDrawing();
    }
    CloseWindow();
    fclose(fp);
}
