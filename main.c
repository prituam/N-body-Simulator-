#include "raylib.h"
#include <math.h>
#define RAYGUI_IMPLEMENTATION
#include "raygui.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
const int G = 10000;

#define MAX_TRAIL 3000 // Maximum number of trail points
int mainmenu = 0;
typedef enum
{
    PLAY,
    PAUSE
} state;

state gamestate = PAUSE;

typedef struct
{
    double x, y;
    double vx, vy;
} body;

typedef struct
{
    double x, y;
} Force;

typedef struct
{
    float x, y, r, a;
} star;
int trail = 1;
Force computeForce(body a, body b)
{
    double dx = b.x - a.x;
    double dy = b.y - a.y;
    double r = sqrt(dx * dx + dy * dy); // distance between the centres of the planets
    double mag = G / (r * r * r);       /// magnitude of the force (vector not included)
    Force f;
    f.x = mag * dx;
    f.y = mag * dy;
    return f;
}

void collide(body *a, body *b)
{
    double dx = b->x - a->x;
    double dy = b->y - a->y;
    double dist2 = dx * dx + dy * dy;
    if (dist2 == 0)
        dist2 = 0.01;

    double dvx = a->vx - b->vx;
    double dvy = a->vy - b->vy;

    double dot = (dvx * dx + dvy * dy) / dist2;

    double vx_a = a->vx - dot * dx;
    double vy_a = a->vy - dot * dy;
    double vx_b = b->vx + dot * dx;
    double vy_b = b->vy + dot * dy;

    a->vx = vx_a;
    a->vy = vy_a;
    b->vx = vx_b;
    b->vy = vy_b;

    // Push apart to avoid sticking
    double dist = sqrt(dist2);
    double overlap = 20 - dist;
    if (overlap > 0)
    {
        double nx = dx / dist;
        double ny = dy / dist;
        a->x -= nx * overlap / 2;
        a->y -= ny * overlap / 2;
        b->x += nx * overlap / 2;
        b->y += ny * overlap / 2;
    }
}
int t = 0; // in order to change the alpha value of the stars

int main()
{
    Color black = {0, 0, 0, 200};
    body b[4];
    star s[200];
    InitWindow(1280, 720, "FOUR STAR SYSTEM WITH TRAILS");
    SetTargetFPS(170);

    // calculating randoom positions for the 200 stars
    for (int i = 0; i < 200; i++)
    {
        s[i].x = GetRandomValue(0, 1280);
        s[i].y = GetRandomValue(0, 720);
        s[i].a = 0;
    }
    // Initial positions
    b[0].x = 540;
    b[0].y = 260;
    b[1].x = 540;
    b[1].y = 460;
    b[2].x = 740;
    b[2].y = 260;
    b[3].x = 740;
    b[3].y = 46;

    // setting initial velocities of the planet 0
    for (int i = 0; i < 4; i++)
    {
        b[i].vx = b[i].vy = 0;
    }

    body rb[4];

    Force force[4] ;
    for(int i=0;i<4;i++)
    {
        force[i].x=0;
        force[i].y=0;
    }

    // Trails array
    Vector2 trails[4][MAX_TRAIL];
    int trailIndex = 0;

    while (!WindowShouldClose())
    {

        BeginDrawing();
        if (mainmenu)
        {
            double dt = GetFrameTime();
            // draw random stars at start of the window
            for (int i = 0; i < 200; i++)

            {
                if (s[i].a >= 255)
                    t = 1;
                if (s[i].a <= 0)
                    t = 0;
                if (t)
                    s[i].a -= 20;
                else
                    s[i].a += 20;
            }
            for (int i = 0; i < 200; i++)
            {
                Color c = {240, 240, 255, s[i].a};
                DrawPixel(s[i].x, s[i].y, c);
            }
            // Reset forces
            if (!gamestate)
            {
                
                for (int i = 0; i < 4; i++)
                {
                    force[i].x = 0;
                    force[i].y = 0;
                }

                // Compute forces
                for (int i = 0; i < 4; i++)
                {
                    for (int j = i + 1; j < 4; j++)
                    {
                        Force f = computeForce(b[i], b[j]);
                        force[i].x += f.x;
                        force[i].y += f.y;
                        force[j].x -= f.x;
                        force[j].y -= f.y;
                    }
                }

                // Update velocities
                for (int i = 0; i < 4; i++)
                {
                    b[i].vx += force[i].x * dt;
                    b[i].vy += force[i].y * dt;
                }

                // collision detection and aftermath
                for (int i = 0; i < 4; i++)
                {
                    for (int j = i + 1; j < 4; j++)
                    {
                        double dx = b[i].x - b[j].x;
                        double dy = b[i].y - b[j].y;
                        double r = sqrt(dx * dx + dy * dy);
                        if (r < 20)
                        {
                            collide(&b[i], &b[j]);
                        }
                    }
                }
                // Update positions
                for (int i = 0; i < 4; i++)
                {
                    b[i].x += b[i].vx * dt;
                    b[i].y += b[i].vy * dt;
                }
            }
            // Bounce from walls
            for (int i = 0; i < 4; i++)
            {
                if (b[i].x > 1270 || b[i].x < 10)
                    b[i].vx = -b[i].vx;
                if (b[i].y > 710 || b[i].y < 10)
                    b[i].vy = -b[i].vy;
            }

            // Store trails
            for (int i = 0; i < 4; i++)
            {
                trails[i][trailIndex] = (Vector2){(float)b[i].x, (float)b[i].y};
            }
            trailIndex = (trailIndex + 1) % MAX_TRAIL;

            // Drawing
            
            // DrawRectangle for glowing trail and fading trail effect

            DrawRectangle(0, 0, 1280, 720, Fade(BLACK, 0.2f));

            if (GuiButton((Rectangle){20, 20, 80, 30}, "Play"))
            {
                gamestate = PLAY;
            }

            if (GuiButton((Rectangle){110, 20, 80, 30}, "Pause"))
            {
                gamestate = PAUSE;
            }
            if (GuiButton((Rectangle){200, 20, 80, 30}, "RESET"))
            {
                gamestate = PAUSE;
                b[0].x = 540;
                b[0].y = 260;
                b[1].x = 540;
                b[1].y = 460;
                b[2].x = 740;
                b[2].y = 260;
                b[3].x = 740;
                b[3].y = 460;
                for (int i = 0; i < 4; i++)
                {
                    b[i].vx = b[i].vy = 0;
                }
                DrawRectangle(0, 0, 1280, 720, Fade(BLACK, 255));
                trailIndex = 0;
            }
            if (GuiButton((Rectangle){20, 60, 80, 30}, "Show trail"))
            {
                trail = 1;
            }
            if (GuiButton((Rectangle){110, 60, 80, 30}, "Hide Trail"))
            {
                trail = 0;
                trailIndex = 0;
            }

            // Draw trails
            if (trail)
            {
                for (int i = 0; i < 4; i++)
                {
                    Color c = (i == 0 ? WHITE : (i == 1 ? RED : (i == 2 ? BLUE : GREEN)));
                    for (int t = 0; t < MAX_TRAIL; t++)
                    {
                        DrawPixelV(trails[i][t], c);
                    }
                }
            }

            // Draw planets

            DrawCircle(b[0].x, b[0].y, 10, WHITE);
            DrawCircle(b[1].x, b[1].y, 10, RED);
            DrawCircle(b[2].x, b[2].y, 10, BLUE);
            DrawCircle(b[3].x, b[3].y, 10, GREEN);

            // twinkling effect on stars
            
        }

        else{
            BeginDrawing();
            DrawText("WELCOME",600,300,30,WHITE);
            if(GuiButton((Rectangle){600,360,80,30},"START"))
            mainmenu=1;
            EndDrawing();
        }

            EndDrawing();
    }
    CloseWindow();
    return 0;
}