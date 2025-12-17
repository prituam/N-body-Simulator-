#include <stdio.h>
#define RAYGUI_IMPLEMENTATION
#include "raygui.h"
#include <math.h>
int mainmenu = 1;
int integrator;

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
    double r = sqrt(dx * dx + dy * dy)+1e-6;
    double mag = (G * a.m * b.m) / (r * r * r);
    Force f;
    f.x = mag * dx;
    f.y = mag * dy;
    return f;
}

void updatePhysics(int i, int integrator, Force f, double dt, Body p[])
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
    case 3: // rk2
        // a.vx+=f.x*dt/2;
        Force ftemp = {0, 0};
        double ax = f.x / p[i].m;
        double ay = f.y / p[i].m;
        double xtemp = p[i].x + p[i].vx * dt / 2;
        double vxtemp = p[i].vx + ax * dt / 2;
        double ytemp = p[i].y + p[i].vy * dt / 2;
        double vytemp = p[i].vy + ay * dt / 2;
        Body tempbody = {xtemp, ytemp, vxtemp, vytemp, p[i].m, 10};

        for (int j = 0; j < 9; j++)
        {
            if (j == i)
                continue;
            else
            {
                ftemp.x += (computeForce(tempbody, p[j])).x;
                ftemp.y += (computeForce(tempbody, p[j])).y;
            }
        }
        double axmid = ftemp.x / p[i].m;
        double aymid = ftemp.y / p[i].m;
        p[i].vx += axmid * dt;
        p[i].vy += aymid * dt;
        p[i].x += (p[i].vx + vxtemp) * 0.5 * dt; // optional
        p[i].y += (p[i].vy + vytemp) * 0.5 * dt;
        break;
    case 4: // RK4
    {
        // --------------------------------------------------
        // STEP 0: original state (anchor)
        // --------------------------------------------------
        double x0 = p[i].x;
        double y0 = p[i].y;
        double vx0 = p[i].vx;
        double vy0 = p[i].vy;
        double m = p[i].m;

        // --------------------------------------------------
        // k1: slopes at original state
        // --------------------------------------------------
        double k1_x = vx0;
        double k1_y = vy0;
        double k1_vx = f.x / m;
        double k1_vy = f.y / m;

        // --------------------------------------------------
        // k2: slopes at half-step using k1
        // --------------------------------------------------
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

        // --------------------------------------------------
        // k3: slopes at half-step using k2
        // --------------------------------------------------
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

        // --------------------------------------------------
        // k4: slopes at full-step using k3
        // --------------------------------------------------
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

        // --------------------------------------------------
        // FINAL UPDATE (single real step)
        // --------------------------------------------------
        p[i].x += (dt / 6) * (k1_x + 2 * k2_x + 2 * k3_x + k4_x);
        p[i].y += (dt / 6) * (k1_y + 2 * k2_y + 2 * k3_y + k4_y);
        p[i].vx += (dt / 6) * (k1_vx + 2 * k2_vx + 2 * k3_vx + k4_vx);
        p[i].vy += (dt / 6) * (k1_vy + 2 * k2_vy + 2 * k3_vy + k4_vy);

        break;
    }

    case 5: // leepfrog
        break;

    default:
        break;
    }
}

int main()
{
    InitWindow(1280, 720, "Solar System");
    SetTargetFPS(190);
    Body planet[9];
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
    int sun_index = 8; // Sun's index

    for (int i = 0; i < 8; i++)
    {
        double dx = planet[i].x - planet[sun_index].x;
        double r = fabs(dx); // horizontal distance

        // velocity magnitude for circular orbit
        double v = sqrt(G * planet[sun_index].m / r);

        // set vertical velocity
        if (dx < 0)
            planet[i].vy = v; // left of Sun, orbit clockwise
        else
            planet[i].vy = -v; // right of Sun, orbit counterclockwise

        // horizontal velocity = 0
        planet[i].vx = 0;
    }

    Force force[9];
    double dt;
    while (!WindowShouldClose())
    {
        dt = GetFrameTime();
        BeginDrawing();
        ClearBackground(RAYWHITE);
        if (mainmenu)
        {
            DrawText("Please Select Numerical Integrator", 500, 220, 19, WHITE);
            if (GuiButton((Rectangle){610, 250, 80, 30}, "Euler"))
            {
                mainmenu = 0;
                integrator = 1;
            }
            else if (GuiButton((Rectangle){590, 285, 120, 30}, "Simplectic Euler"))
            {
                mainmenu = 0;
                integrator = 2;
            }
            else if (GuiButton((Rectangle){610, 320, 80, 30}, "RK2"))
            {
                mainmenu = 0;
                integrator = 3;
            }
            else if (GuiButton((Rectangle){610, 355, 80, 30}, "RK4"))
            {
                mainmenu = 0;
                integrator = 4;
            }
            else if (GuiButton((Rectangle){610, 390, 80, 30}, "Leepfrog"))
            {
                mainmenu = 0;
                integrator = 5;
            }
        }
        // main physics and drawing starts from here
        else
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
            for (int i = 0; i < 9; i++)
            {
                updatePhysics(i, integrator, force[i], dt, planet);
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
