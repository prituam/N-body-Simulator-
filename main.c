#include <stdio.h>
#define RAYGUI_IMPLEMENTATION
#include "raygui.h"
#include <math.h>
int mainmenu = 1;
int integrator;
int play = 1;
const double G = 0.2;

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
    double r = sqrt(dx * dx + dy * dy) + 1e-6;
    double mag = (G * a.m * b.m) / (r * r * r);
    Force f;
    f.x = mag * dx;
    f.y = mag * dy;
    return f;
}
void applyRK2(double dt, Body p[], Force f[])
{
    Force fmid[9];
    Body temp[9];
    double vx_mid[9], vy_mid[9];


    for (int i = 0; i < 9; i++)
        fmid[i] = (Force){0, 0};


    for (int i = 0; i < 9; i++)
    {
        double ax = f[i].x / p[i].m;
        double ay = f[i].y / p[i].m;

        temp[i].x  = p[i].x  + p[i].vx * dt * 0.5;
        temp[i].y  = p[i].y  + p[i].vy * dt * 0.5;
        temp[i].vx = p[i].vx + ax * dt * 0.5;
        temp[i].vy = p[i].vy + ay * dt * 0.5;
        temp[i].m  = p[i].m;

        vx_mid[i] = temp[i].vx;
        vy_mid[i] = temp[i].vy;
    }

 
    for (int i = 0; i < 9; i++)
        for (int j = i + 1; j < 9; j++)
        {
            Force fij = computeForce(temp[i], temp[j]);
            fmid[i].x += fij.x;
            fmid[i].y += fij.y;
            fmid[j].x -= fij.x;
            fmid[j].y -= fij.y;
        }


    for (int i = 0; i < 9; i++)
    {
        p[i].vx += (fmid[i].x / p[i].m) * dt;
        p[i].vy += (fmid[i].y / p[i].m) * dt;
        p[i].x  += vx_mid[i] * dt;
        p[i].y  += vy_mid[i] * dt;
    }
}

void updatePhysics(int i, int integrator, Force f, double dt, Body p[], Force fo[])
{
    switch (integrator)
    {
    case 1: // euler
        p[i].x += p[i].vx * dt;
        p[i].y += p[i].vy * dt;
        p[i].vx += (f.x / p[i].m) * dt;
        p[i].vy += (f.y / p[i].m) * dt;
        break;
    case 2: // simplectic euler------------------------
        p[i].vx += (f.x / p[i].m) * dt;
        p[i].vy += (f.y / p[i].m) * dt;
        p[i].x += p[i].vx * dt;
        p[i].y += p[i].vy * dt;
        break;
    case 4: // RK4
    {

        double x0 = p[i].x;
        double y0 = p[i].y;
        double vx0 = p[i].vx;
        double vy0 = p[i].vy;
        double m = p[i].m;

        double k1_x = vx0;
        double k1_y = vy0;
        double k1_vx = f.x / m;
        double k1_vy = f.y / m;

        Body s2;
        s2.x = x0 + k1_x * dt / 2;
        s2.y = y0 + k1_y * dt / 2;
        s2.vx = vx0 + k1_vx * dt / 2;
        s2.vy = vy0 + k1_vy * dt / 2;
        s2.m = m;

        Force f2 = {0, 0};
        for (int j = 0; j < 9; j++)
        {
            if (j != i)
            {
                f2.x += (computeForce(s2, p[j])).x;
                f2.y += (computeForce(s2, p[j])).y;
            }
        }

        double k2_x = s2.vx;
        double k2_y = s2.vy;
        double k2_vx = f2.x / m;
        double k2_vy = f2.y / m;

        Body s3;
        s3.x = x0 + k2_x * dt / 2;
        s3.y = y0 + k2_y * dt / 2;
        s3.vx = vx0 + k2_vx * dt / 2;
        s3.vy = vy0 + k2_vy * dt / 2;
        s3.m = m;

        Force f3 = {0, 0};
        for (int j = 0; j < 9; j++)
        {
            if (j != i)
            {
                f3.x += (computeForce(s3, p[j])).x;
                f3.y += (computeForce(s3, p[j])).y;
            }
        }

        double k3_x = s3.vx;
        double k3_y = s3.vy;
        double k3_vx = f3.x / m;
        double k3_vy = f3.y / m;

        Body s4;
        s4.x = x0 + k3_x * dt;
        s4.y = y0 + k3_y * dt;
        s4.vx = vx0 + k3_vx * dt;
        s4.vy = vy0 + k3_vy * dt;
        s4.m = m;

        Force f4 = {0, 0};
        for (int j = 0; j < 9; j++)
        {
            if (j != i)
            {
                f4.x += (computeForce(s4, p[j])).x;
                f4.y += (computeForce(s4, p[j])).y;
            }
        }

        double k4_x = s4.vx;
        double k4_y = s4.vy;
        double k4_vx = f4.x / m;
        double k4_vy = f4.y / m;

        p[i].x += (dt / 6) * (k1_x + 2 * k2_x + 2 * k3_x + k4_x);
        p[i].y += (dt / 6) * (k1_y + 2 * k2_y + 2 * k3_y + k4_y);
        p[i].vx += (dt / 6) * (k1_vx + 2 * k2_vx + 2 * k3_vx + k4_vx);
        p[i].vy += (dt / 6) * (k1_vy + 2 * k2_vy + 2 * k3_vy + k4_vy);

        break;
    }

    default:
        break;
    }
}
void applyLeapfrog(double dt, Body p[], Force f[])
{

    for (int i = 0; i < 9; i++)
    {
        p[i].vx += (f[i].x / p[i].m) * dt * 0.5;
        p[i].vy += (f[i].y / p[i].m) * dt * 0.5;

        p[i].x += p[i].vx * dt;
        p[i].y += p[i].vy * dt;
    }

    for (int i = 0; i < 9; i++)
        f[i] = (Force){0, 0};

    for (int i = 0; i < 9; i++)
    {
        for (int j = i + 1; j < 9; j++)
        {
            Force fij = computeForce(p[i], p[j]);
            f[i].x += fij.x;
            f[i].y += fij.y;
            f[j].x -= fij.x;
            f[j].y -= fij.y;
        }
    }

    for (int i = 0; i < 9; i++)
    {
        p[i].vx += (f[i].x / p[i].m) * dt * 0.5;
        p[i].vy += (f[i].y / p[i].m) * dt * 0.5;
    }
}

void initialize(Body planet[])
{
    planet[0].x = 648;  // mercury
    planet[1].x = 654;  // venus
    planet[2].x = 660;  // earth
    planet[3].x = 670;  // mars
    planet[4].x = 744;  // jupiter
    planet[5].x = 832;  // saturn
    planet[6].x = 1024; // uranus
    planet[7].x = 1241; // neptune
    planet[8].x = 640;  // sun

    planet[0].m = 0.0000165; // mercury
    planet[1].m = 0.000245;  // venus
    planet[2].m = 0.0003;    // earth
    planet[3].m = 0.000032;  // mars
    planet[4].m = 0.0954;    // jupiter
    planet[5].m = 0.0544;    // saturn
    planet[6].m = 0.0081;    // uranus
    planet[7].m = 0.0102;    // neptune
    planet[8].m = 100;       // sun

    for (int i = 0; i < 9; i++)
    {
        planet[i].y = 360;
        planet[i].vx = 0;
    }
    // scaled gravitational constant
    int sun_index = 8;

    for (int i = 0; i < 8; i++)
    {
        double dx = planet[i].x - planet[sun_index].x;
        double r = fabs(dx);

        double v = sqrt(G * planet[sun_index].m / r);

        if (dx < 0)
            planet[i].vy = v;
        else
            planet[i].vy = -v;

        // horizontal velocity = 0
        planet[i].vx = 0;
    }
}

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
    while (!WindowShouldClose())
    {
        dt = 0.002;
        BeginDrawing();
        ClearBackground(RAYWHITE);
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
        // main physics and drawing starts from here
        else
        {
            if (GuiButton((Rectangle){13, 10, 175, 30}, "Change Integrator"))
            {
                mainmenu = 1;
                PlaySound(clickSound);
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
                    PlaySound(clickSound);
                    play = 1;
                }
            }
            else
            {
                if (GuiButton((Rectangle){13, 80, 175, 30}, "Pause Simulation"))
                {
                    PlaySound(clickSound);
                    play = 0;
                }
            }
            if (play)
            {
                for (int i = 0; i < 9; i++)
                {
                    force[i].x = 0;
                    force[i].y = 0;
                }
                for (int i = 0; i < 9; i++)
                {
                    for (int j = i + 1; j < 9; j++)
                    {
                        Force f = computeForce(planet[i], planet[j]);
                        force[i].x += f.x;
                        force[i].y += f.y;
                        force[j].x -= f.x;
                        force[j].y -= f.y;
                    }
                }
                if (integrator == 3)
                {
                    applyRK2(dt, planet, force);
                }
                else if (integrator == 5)
                {
                    applyLeapfrog(dt, planet, force);
                }
                else
                {
                    for (int i = 0; i < 9; i++)
                    {
                        updatePhysics(i, integrator, force[i], dt, planet, force);
                    }
                }
            }
            for (int i = 0; i < 9; i++)
            {
                DrawCircle(planet[i].x, planet[i].y, 2, GRAY);
            }
        }
        EndDrawing();
    }
    CloseWindow();
}
