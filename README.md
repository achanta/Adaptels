# Adaptels

This repository offers the code for the Adaptels algorithm with a Python interface. A demo file is provided, which should be easy to use. Both versions can produce superpixels for grey, color, as well as images with any other number of channels. If you use the code, please cite the following publication:

"Scale-Adaptive Superpixels", R. Achanta and P Marquez-Neila and P. Fua and S. Süsstrunk, Color and Imaging Conference (CIC26), 2018

```
@article{adaptels_2018,
      title = {Scale-Adaptive Superpixels},
      author = {Achanta, Radhakrishna and Marquez Neila, Pablo and Fua,  Pascal and Süsstrunk, Sabine},
      publisher = {Society for Imaging Science and Technology},
      journal = {26th Color and Imaging Conference Final Program and  Proceedings},
      pages = {1-6(6)},
      year = {2018},
}
```

The Adaptels algorithm is different from conventional superpixel algorithms in that it generates compact segments that adapt to the local image scale. In regions of low information it generates large superpixels and vice versa. It simple to use since it takes just one parameter (which can be set to 60 for sRGB images). It is also computationally very efficient. An example segmentation of an image is shown below:

<p float="center">
  <img src="https://github.com/achanta/Adaptels/blob/master/bee.png" width="400" />
  <img src="https://github.com/achanta/Adaptels/blob/master/bee_adaptels.png" width="400" /> 
</p>
