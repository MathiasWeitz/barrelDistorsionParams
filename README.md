barrelDistorsionParams
======================

calculate barrel distorsion parameter for imagemagick

What it does
------------

Imagemagick provides a barrel distorsion with 4 parameters. While one of these parameters is for linear zooming, the other three are for a polynom distorsion
http://www.imagemagick.org/Usage/distorts/#barrel
This script calculate the 3 parameters by a numeric approximation

Motivation
----------

I have got a Actioncam (Contour1080) with a nasty distorsion. I have tested several programs for detorsion and found them all wanting.
The main problem is, that i hardly get a straight line by pure visual estimation. Another problem is the lavish cropping.
Imagemagick offers a very accurate detorsion, with a good control of cropping, but yoa have to know the parameters of your cam. And if you want to change the cropping, you even have different parameters. Especially for me, i noticed, that i get much more of the original image if i allow a little from the outside image to appear on the canvas, also i can crop the top and the bottom of the image and therefore change the ratio of the film.
All this requires that i determine the correct parameters.

The only program (i found) to determine the parameter is Hugin, but it is complicated to operate, take a lot of time, quite often fails with a result, and even if it delievers a result, it is sometimes not accurate. Also, i can not determine the parameter by a given zoom factor.

What you need to do
-------------------

The whole program runs on Python, there is no GUI, you have to write the data right into the code.
So a little knowledge of Python is helping.

In the program you will find a array with array of coordinates. These coordinates are on the image what is in real world a straight line - but on the image distorted. I have determined this coordinates with gimp and wrote them down. A straight line should at least contain three coordinates, more are of course better. You can write down straight lines from different images.
You also have to enter the size of the canvas, the program uses is for normalisation of the coordinates.

You have a function @calculateDistorsion@ which two parameters, the first one are the coordinates, the second is the zoom (number goes up, you zoom out).

What you get
------------

When running the program from the console you will see the program in several steps approximate the parameters.
You get the number of the steps, then an estimate of the value of correction (it should approx zero) and a value of error (it should be small)

At the end you get the four values a,b,c,d where d is the value you have put in.

The program may fail due to the instability to the numeric approximation, in this case the value of correction remains high.


