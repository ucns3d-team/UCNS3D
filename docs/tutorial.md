
# 2D Solid Body Rotation 


table th {
         background-color: #FF0000;
        /* sets table header cell background colour */
    }
  
| Written by       | for Version | Revised by | Revision Date |for Version |
| :---------------: | :------: | :---: | :------: | :------: |
| Python Hat        |   True   | 23.99 | True   | 23.99 |


## Summary

![Alt text](SBR2.gif)


The solid body rotation test of Leveque \cite{Leveque1996627} is employed to investigate the performance of the WENO, CWENO and CWENOZ schemes in terms of their non-oscillatory properties on the correct and accurate revolution of a body

## Aim

## Resources

## Tutorial

## Background

## Problem Setup

## Mesh Description


## Configuration File Options


## Running UCNS3D

## Results

## Description
The solid body rotation test of Leveque \cite{Leveque1996627} is employed to investigate the performance of the WENO, CWENO and CWENOZ schemes in terms of their non-oscillatory properties on the correct and accurate revolution of a body. The continuity equation is considered as follows:


$\frac{\partial{U}}{\partial t}+\nabla \cdot({\mathbf{v}U})=0$

$\mathbf{v}(x,y)=(0.5-y,x-0.5)$


In this particular test three bodies are considered namely a smooth centered at $(x_0=0.25,y_0=0.5)$, a sharp cone centered at $(x_0=0.5,y_0=0.25)$ and a slotted cylinder centered at $(x_0=0.5,y_0=0.75)$ described by the following functions respectively:


$\mathbf{f}(x,y)=\frac{1+cos(\pi r(x,y))}{4}$


$\mathbf{f}(x,y)=1- r(x,y)$


with $r_0=0.15$. The rest of the domain the solution is initialised with zero, and after one full revolution $t_{f}=2\pi$ the exact solution coincides with the initial solution. A triangular unstructured mesh is used as shown in  with $64$ edges per side of the computational domain. The WENO, CWENO and CWENOZ schemes ranging from 3rd- to 6th-order of spatial accuracy are employed, with their selected weights from the previous study. The obtained results following one revolution are illustrated in Fig. .

## Step 1


## Step 2


## Step 3


## Step 4

Bullet list:
* Start a line with an asterisk
* Food
  * Fruits
    * Oranges
    * Apples

Numbered list:
1. One
2. Two
3. Three

If you have inline code blocks, you can wrap them in backticks: `var example = true`.

If you've got a longer block of code, you can use ```:

```
if (isAwesome){
  return true
}
```



