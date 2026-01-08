#include <stdio.h>
#include <math.h>
#define RAYGUI_IMPLEMENTATION
#include "raygui.h"
#define GRAPH_N 600
int d = 0;
double energyarr[GRAPH_N];
int idx = 0;
double scale = 1e11;
const int x0 = 900;
const int yo = 220;
const double scaleE = 1e-31; // adjustable
int planetselected = -1;
int mainmenu = 1;
int integrator;
int play = 1;
const double G = 20;
const double energyFactor = 6.6e33;
int leapfroginit = 0;
int energyinit = 0;
int showVelVector=1;
typedef struct
{
    double x, y, vx, vy, m, r;
    Color c;
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
double KE(Body p[], Force f[], int integrator, double dt)
{
    double ke = 0.0;

    for (int i = 0; i < 9; i++)
    {
        double vx, vy;

        vx = p[i].vx;
        vy = p[i].vy;

        ke += 0.5 * p[i].m * (vx * vx + vy * vy);
    }
    return ke;
}

double PE(Body p[])
{
    double pe = 0.0;
    for (int i = 0; i < 9; i++)
        for (int j = i + 1; j < 9; j++)
        {
            double dx = p[j].x - p[i].x;
            double dy = p[j].y - p[i].y;
            double r = sqrt(dx * dx + dy * dy + 1e-6);
            pe -= G * p[i].m * p[j].m / r;
        }
    return pe;
}

double calcTE(Body p[], Force force[], double dt)
{
    return KE(p, force, integrator, dt) + PE(p);
}

void applyRK2(double dt, Body p[], Force f[])
{
    DrawText("RK2", x0 + 120, yo - 150, 12, WHITE);
    scale = 1e13;
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
        p[i].x += temp[i].vx * dt;
        p[i].y += temp[i].vy * dt;

        // Velocities updated using **mid-step accelerations**
        p[i].vx += ax_mid * dt;
        p[i].vy += ay_mid * dt;
    }
}
void applyRK4(double dt, Body p[], Force f[])
{
    DrawText("RK4", x0 + 120, yo - 150, 12, WHITE);
    scale = 1e15;
    // double x0 = p[i].x, y0 = p[i].y;
    // double vx0 = p[i].vx, vy0 = p[i].vy;
    double k1_x[9];
    double k1_vx[9];
    double k1_y[9];
    double k1_vy[9];
    Force f1[9] = {0};
    for (int i = 0; i < 9; i++)
        for (int j = i + 1; j < 9; j++)
        {
            Force fij = computeForce(p[i], p[j]);
            f1[i].x += fij.x;
            f1[i].y += fij.y;
            f1[j].x -= fij.x;
            f1[j].y -= fij.y;
        }

    for (int i = 0; i < 9; i++)
    {
        k1_x[i] = p[i].vx;
        k1_vx[i] = f1[i].x / p[i].m;
        k1_y[i] = p[i].vy;
        k1_vy[i] = f1[i].y / p[i].m;
    }
    Body s2[9];
    for (int i = 0; i < 9; i++)
    {
        s2[i] = (Body){p[i].x + k1_x[i] * dt / 2, p[i].y + k1_y[i] * dt / 2, p[i].vx + k1_vx[i] * dt / 2, p[i].vy + k1_vy[i] * dt / 2, p[i].m, 0};
    }

    Force f2[9];
    for (int i = 0; i < 9; i++)
    {
        f2[i] = (Force){0, 0};
    }

    for (int i = 0; i < 9; i++)
        for (int j = i + 1; j < 9; j++)
        {
            Force fx = computeForce(s2[i], s2[j]);
            f2[i].x += fx.x;
            f2[i].y += fx.y;
            f2[j].x -= fx.x;
            f2[j].y -= fx.y;
        }

    double k2_x[9];
    double k2_vx[9];
    double k2_y[9];
    double k2_vy[9];
    for (int i = 0; i < 9; i++)
    {
        k2_x[i] = s2[i].vx;
        k2_vx[i] = f2[i].x / s2[i].m;
        k2_y[i] = s2[i].vy;
        k2_vy[i] = f2[i].y / s2[i].m;
    }
    Body s3[9];
    for (int i = 0; i < 9; i++)
    {
        s3[i] = (Body){p[i].x + k2_x[i] * dt / 2, p[i].y + k2_y[i] * dt / 2, p[i].vx + k2_vx[i] * dt / 2, p[i].vy + k2_vy[i] * dt / 2, p[i].m, 0};
    }
    Force f3[9];
    for (int i = 0; i < 9; i++)
    {
        f3[i] = (Force){0, 0};
    }

    for (int i = 0; i < 9; i++)
        for (int j = i + 1; j < 9; j++)
        {
            Force fx = computeForce(s3[i], s3[j]);
            f3[i].x += fx.x;
            f3[i].y += fx.y;
            f3[j].x -= fx.x;
            f3[j].y -= fx.y;
        }

    double k3_x[9];
    double k3_vx[9];
    double k3_y[9];
    double k3_vy[9];
    for (int i = 0; i < 9; i++)
    {
        k3_x[i] = s3[i].vx;
        k3_vx[i] = f3[i].x / s3[i].m;
        k3_y[i] = s3[i].vy;
        k3_vy[i] = f3[i].y / s3[i].m;
    }
    Body s4[9];
    for (int i = 0; i < 9; i++)
    {
        s4[i] = (Body){p[i].x + k3_x[i] * dt, p[i].y + k3_y[i] * dt, p[i].vx + k3_vx[i] * dt, p[i].vy + k3_vy[i] * dt, p[i].m, 0};
    }
    Force f4[9];
    for (int i = 0; i < 9; i++)
    {
        f4[i] = (Force){0, 0};
    }

    for (int i = 0; i < 9; i++)
        for (int j = i + 1; j < 9; j++)
        {
            Force fx = computeForce(s4[i], s4[j]);
            f4[i].x += fx.x;
            f4[i].y += fx.y;
            f4[j].x -= fx.x;
            f4[j].y -= fx.y;
        }

    double k4_x[9];
    double k4_vx[9];
    double k4_y[9];
    double k4_vy[9];
    for (int i = 0; i < 9; i++)
    {
        k4_x[i] = s4[i].vx;
        k4_vx[i] = f4[i].x / s4[i].m;
        k4_y[i] = s4[i].vy;
        k4_vy[i] = f4[i].y / s4[i].m;
    }
    for (int i = 0; i < 9; i++)
    {
        p[i].x += dt / 6 * (k1_x[i] + 2 * k2_x[i] + 2 * k3_x[i] + k4_x[i]);
        p[i].y += dt / 6 * (k1_y[i] + 2 * k2_y[i] + 2 * k3_y[i] + k4_y[i]);
        p[i].vx += dt / 6 * (k1_vx[i] + 2 * k2_vx[i] + 2 * k3_vx[i] + k4_vx[i]);

        p[i].vy += dt / 6 * (k1_vy[i] + 2 * k2_vy[i] + 2 * k3_vy[i] + k4_vy[i]);
    }
}
void updatePhysics(int i, int integrator, Force f, double dt, Body p[], Force fo[])
{
    switch (integrator)
    {
    case 1: // Euler
        DrawText("Euler", x0 + 120, yo - 150, 12, WHITE);
        scale = 1e7;
        p[i].x += p[i].vx * dt;
        p[i].y += p[i].vy * dt;
        p[i].vx += (f.x / p[i].m) * dt;
        p[i].vy += (f.y / p[i].m) * dt;

        break;
    case 2: // Symplectic Euler
        scale = 1e9;
        DrawText("Simplectic Euler", x0 + 120, yo - 150, 12, WHITE);
        p[i].vx += (f.x / p[i].m) * dt;
        p[i].vy += (f.y / p[i].m) * dt;
        p[i].x += p[i].vx * dt;
        p[i].y += p[i].vy * dt;

        break;
    default:
        break;
    }
}

void applyLeapfrog(double dt, Body p[], Force f[])
{
    DrawText("Leapfrog(short x axis)", x0 + 120, yo - 150, 12, WHITE);
    scale = 1e12;
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

    planet[0].x = 600; // Mercury

    planet[0].vx = 0;
    planet[0].vy = 0;
    planet[0].m = 0.0000165;
    planet[0].r = 4;

    planet[1].x = 630; // Venus

    planet[1].vx = 0;
    planet[1].vy = 0;
    planet[1].m = 0.000245;
    planet[1].r = 3;

    planet[2].x = 660; // Earth

    planet[2].vx = 0;
    planet[2].vy = 0;
    planet[2].m = 0.0003;
    planet[2].r = 5.5;

    planet[3].x = 700; // Mars

    planet[3].vx = 0;
    planet[3].vy = 0;
    planet[3].m = 0.000032;
    planet[3].r = 6;

    planet[4].x = 760; // Jupiter

    planet[4].vx = 0;
    planet[4].vy = 0;
    planet[4].m = 0.0954;
    planet[4].r = 9;

    planet[5].x = 860; // Saturn

    planet[5].vx = 0;
    planet[5].vy = 0;
    planet[5].m = 0.0544;
    planet[5].r = 8;

    planet[6].x = 1020; // Uranus

    planet[6].vx = 0;
    planet[6].vy = 0;
    planet[6].m = 0.0081;
    planet[6].r = 6;

    planet[7].x = 1200; // Neptune

    planet[7].vx = 0;
    planet[7].vy = 0;
    planet[7].m = 0.0102;
    planet[7].r = 6;

    planet[8].x = 640;

    planet[8].vx = 0;
    planet[8].vy = 0;
    planet[8].m = 100;
    planet[8].r = 6;
    planet[8].y = 360;

    planet[0].c = (Color){255, 255, 204, 255}; // mercury
    planet[1].c = (Color){255, 128, 0, 255};   // venus
    planet[2].c = (Color){51, 255, 255, 255};  // earth
    planet[3].c = (Color){255, 204, 153, 255}; // mars
    planet[4].c = (Color){255, 102, 102, 255}; // jupiter
    planet[5].c = (Color){255, 229, 204, 255}; // saturn
    planet[6].c = (Color){102, 178, 255, 255}; // uranus
    planet[7].c = (Color){102, 102, 255, 255}; // neptune
    planet[8].c = (Color){255, 255, 153, 255}; // sun

    for (int i = 0; i < 8; i++)
    {
        double dx = planet[i].x - planet[sun_index].x;
        double r = fabs(dx);
        double v = sqrt(G * planet[sun_index].m / r);

        planet[i].vy = ((dx < 0) ? v : -v);
        planet[i].y = 360 + i * 0.5;
    }

    double Px = 0, Py = 0;
    for (int i = 0; i < 8; i++)
    {
        Px += planet[i].m * planet[i].vx;
        Py += planet[i].m * planet[i].vy;
    }
    planet[sun_index].vx = -Px / planet[sun_index].m;
    planet[sun_index].vy = -Py / planet[sun_index].m;
}

void showvector(Body a)
{
    float L = 8.0f;
    float xa = a.x + L * a.vx;
    float ya = a.y + L * a.vy;
    float dx = a.vx;
    float dy = a.vy;
    DrawLine(a.x, a.y, xa, ya, WHITE);
    float h = 1.0f;
    float c = cosf(PI / 6.0f);
    float s = sinf(PI / 6.0f);
    float x1 = xa - h * (dx * c - dy * s);
    float y1 = ya - h * (dx * s + dy * c);

    float x2 = xa - h * (dx * c + dy * s);
    float y2 = ya - h * (-dx * s + dy * c);
    DrawLine(xa, ya, x1, y1, WHITE); // left head
    DrawLine(xa, ya, x2, y2, WHITE);
}
double vecdist(Body a, Vector2 b)
{
    double dx = a.x - b.x;
    double dy = a.y - b.y;
    return sqrt(dx * dx + dy * dy);
}
int main()
{
    int x = 1000;
    double TE0;
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
    FILE *fp = fopen("Simulation_data.csv", "w");
    fprintf(fp, "i,time,x,y,vx,vy,KE,PE,TE\n");
    double t = 0;
    while (!WindowShouldClose())
    {
        x = 1000;
        dt = 0.0002;
        BeginDrawing();
        ClearBackground(BLACK);

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
            DrawText(TextFormat("Energy: %.8e", TE), 30, 140, 16, WHITE);
            if (GuiButton((Rectangle){13, 10, 175, 30}, "Change Integrator"))
            {
                idx = 0;
                mainmenu = 1;
                leapfroginit = 0;
                energyinit = 0;
                PlaySound(clickSound);
                ClearBackground(BLACK);
                initialize(planet);
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
            }
            if (GuiButton((Rectangle){13, 45, 175, 30}, "Reset Simulation"))
            {
                d = 0;
                idx = 0;
                leapfroginit = 0;
                energyinit = 0;

                ClearBackground(BLACK);
                initialize(planet);
                PlaySound(clickSound);
            }

            if (GuiButton((Rectangle){190, 10, 175, 30}, "Rescale Graph"))
            {
                ClearBackground(BLACK);
                TE0 = TE;
                for (int i = 0; i < GRAPH_N - 1; i++)
                {
                    energyarr[i] = 0;
                }
            }
            if (!play)
            {
                if (GuiButton((Rectangle){190, 80, 175, 30}, "Play Simulation"))
                {
                    play = 1;
                    PlaySound(clickSound);
                }
            }
            else
            {
                if (GuiButton((Rectangle){190, 80, 175, 30}, "Pause Simulation"))
                {
                    play = 0;
                    PlaySound(clickSound);
                }
            }
            if (!showVelVector)
            {
                if (GuiButton((Rectangle){13, 80, 175, 30}, "Show Velocity Vector"))
                {
                    showVelVector = 1;
                    PlaySound(clickSound);
                }
            }
            else
            {
                if (GuiButton((Rectangle){13, 80, 175, 30}, "Hide Velocity Vector"))
                {
                    showVelVector = 0;
                    PlaySound(clickSound);
                }
            }
            Vector2 mouse = GetMousePosition();
            int r = 64;
            int ry = 36;
            int minIndex = -1;
            float minDist = 1e30; // very large number

            for (int i = 0; i < 9; i++)
            {
                float d = vecdist(planet[i], mouse);
                if (d < minDist)
                {
                    minDist = d;
                    minIndex = i;
                }
            }
            if (IsMouseButtonPressed(MOUSE_LEFT_BUTTON))
            {
                if (minIndex != -1 && minDist <= planet[minIndex].r)
                {
                    planetselected = 1;
                }
            }
            if(showVelVector)
            if (planetselected)
                showvector(planet[minIndex]);
            for (int i = 0; i < 19; i++)
            {
                DrawLine(0, ry, 1280, ry, (Color){255, 255, 255, 40});
                DrawLine(r, 0, r, 720, (Color){255, 255, 255, 40});
                r += 64;
                ry += 36;
            }
            if (play)
            {
                if (integrator != 5)
                {
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
                }
                if (integrator == 3)
                    applyRK2(dt, planet, force);

                else if (integrator == 5)
                {
                    if (!leapfroginit)
                    {
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
                        leapfroginit = 1;
                    }
                    applyLeapfrog(dt, planet, force);
                }
                else if (integrator == 4)
                {
                    applyRK4(dt, planet, force);
                }
                else
                    for (int i = 0; i < 9; i++)
                        updatePhysics(i, integrator, force[i], dt, planet, force);

                TE = calcTE(planet, force, dt);
                if (!energyinit)
                {
                    TE0 = TE;
                    energyinit = 1;
                }

                if (energyinit)
                {
                    if (idx >= GRAPH_N)
                    {
                        for (int i = 0; i < GRAPH_N - 1; i++)
                        {
                            energyarr[i] = energyarr[i + 1];
                        }
                        energyarr[GRAPH_N - 1] = TE;
                    }
                    if (idx < GRAPH_N)
                    {
                        energyarr[idx] = TE;
                        idx++;
                    }
                }
                for (int i = 0; i < 9; i++)
                {
                    fprintf(fp, "%d,%.6f,%.6f,%.6f,%.6f,%.6f,%.6e,%.6e,%.6e\n",
                            i, t,
                            planet[i].x, planet[i].y,
                            planet[i].vx, planet[i].vy,
                            KE(planet, force, integrator, dt),
                            PE(planet),
                            KE(planet, force, integrator, dt) + PE(planet));
                }

                t += dt;
            }
            for (int i = 0; i < 9; i++)
            {
                DrawCircle(planet[i].x, planet[i].y, planet[i].r, planet[i].c);
            }
            for (int i = 0; i < idx; i++)
            {
                int x = x0 + i / 3;
                double y = yo + scale * (energyarr[i] - TE0) / TE0;
                DrawPixel(x, y, BLUE);
            }
            DrawText("rel. error in TE", x0 - 100, yo - 100, 12, WHITE);
            DrawLine(x0, yo, x0 + 300, yo, WHITE);
            DrawLine(x0, yo - 200, x0, yo, WHITE);
        }

        EndDrawing();
    }
    CloseWindow();
    fclose(fp);
}
