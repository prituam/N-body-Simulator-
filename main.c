#include "raylib.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

const int G = 900;
#define MAX_TRAIL 100 // Maximum number of trail points

typedef struct {
    double x, y;
    double vx, vy;
} body;

typedef struct {
    double x, y;
} Force;

Force computeForce(body a, body b)
{
    double dx = b.x - a.x;
    double dy = b.y - a.y;
    double r = sqrt(dx*dx + dy*dy) + 5;
    double mag = G / (r*r*r);
    Force f;
    f.x = mag * dx;
    f.y = mag * dy;
    return f;    
}

void collide(body *a, body *b)
{
    double dx = b->x - a->x;
    double dy = b->y - a->y;
    double dist2 = dx*dx + dy*dy;
    if(dist2 == 0) dist2 = 0.01;

    double dvx = a->vx - b->vx;
    double dvy = a->vy - b->vy;

    double dot = (dvx*dx + dvy*dy) / dist2;

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
    if(overlap > 0) {
        double nx = dx / dist;
        double ny = dy / dist;
        a->x -= nx * overlap / 2;
        a->y -= ny * overlap / 2;
        b->x += nx * overlap / 2;
        b->y += ny * overlap / 2;
    }
}


int main()
{
    body b[4];
    InitWindow(1280, 720, "FOUR STAR SYSTEM WITH TRAILS");
    SetTargetFPS(150);

    // Initial positions
    b[0].x = 300; b[0].y = 100;
    b[1].x = 200; b[1].y = 200;
    b[2].x = 250; b[2].y = 200;
    b[3].x = 300; b[3].y = 200;

    for(int i = 0; i < 4; i++) {
        b[i].vx = b[i].vy = 0;
    }

    Force force[4] = {0};

    // Trails array
    Vector2 trails[4][MAX_TRAIL];
    int trailIndex = 0;

    while (!WindowShouldClose())
    {
        // Reset forces
        for (int i = 0; i < 4; i++) {
            force[i].x = 0;
            force[i].y = 0;
        }

        // Compute forces
        for (int i = 0; i < 4; i++) {
            for (int j = i + 1; j < 4; j++) {
                Force f = computeForce(b[i], b[j]);
                force[i].x += f.x;
                force[i].y += f.y;
                force[j].x -= f.x;
                force[j].y -= f.y;
            }
        }

        // Update velocities
        for (int i = 0; i < 4; i++) {
            b[i].vx += force[i].x * 0.09;
            b[i].vy += force[i].y * 0.09;
        }

        for (int i = 0; i < 4; i++) {
            for (int j = i + 1; j < 4; j++) {
                double dx=b[i].x-b[j].x;
                double dy=b[i].y-b[j].y;
                double r=sqrt(dx*dx+dy*dy); 
                if(r<20)
                {
                    collide(&b[i],&b[j]);
                }
            }
        }
        // Update positions
        for (int i = 0; i < 4; i++) {
            b[i].x += b[i].vx * 0.09;
            b[i].y += b[i].vy * 0.09;
        }

        // Bounce from walls
        for (int i = 0; i < 4; i++) {
            if (b[i].x > 1270 || b[i].x < 10) b[i].vx = -b[i].vx;
            if (b[i].y > 710 || b[i].y < 10) b[i].vy = -b[i].vy;
        }

        // Store trails
        for (int i = 0; i < 4; i++) {
            trails[i][trailIndex] = (Vector2){(float)b[i].x, (float)b[i].y};
        }
        trailIndex = (trailIndex + 1) % MAX_TRAIL;


        


        // Drawing
        BeginDrawing();
        // DrawRectangle(0, 0, 1280, 720, Fade(BLACK, 0.2f));
        ClearBackground(BLACK);

        // Draw trails
        for (int i = 0; i < 4; i++) {
            Color c = (i == 0 ? WHITE : (i == 1 ? RED : (i == 2 ? BLUE : GREEN)));
            for (int t = 0; t < MAX_TRAIL; t++) {
                DrawPixelV(trails[i][t], c);
            }
        }

        // Draw planets
        DrawCircle(b[0].x, b[0].y, 10, WHITE);
        DrawCircle(b[1].x, b[1].y, 10, RED);
        DrawCircle(b[2].x, b[2].y, 10, BLUE);
        DrawCircle(b[3].x, b[3].y, 10, GREEN);

        EndDrawing();
    }

    CloseWindow();
    return 0;
}
