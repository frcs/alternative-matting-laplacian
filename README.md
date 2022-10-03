# An Alternative Matting Laplacian

This package contains the code described in:

```
[Pitie16] An Alternative Matting Laplacian. F. Pitie, (2016)
          In International Conference on Image Processing (ICIP'16), September.
```

This paper offers an alternative formulation to the Matting Laplacian proposed in [1].

Please cite our publication when using the code.

Send an email to pitief@tcd.ie if you want more information

## Example

Run the demo
```
>>> demo
```

<table style="width:100%">
<tr>
<td><img src="GT04.png"  width="320" ></td>
<td><img src="alpha0-GT04.png"  width="320" ></td>
</tr>
<tr>
<td>input image</td>
<td>input alpha estimate (using [2])</td>
</tr>
<tr>
<td><img src="result-alpha-GT04.png"  width="320" ></td>
<td></td>
</tr>
<tr>
<td>estimated alpha</td>
<td></td>
</tr>
<tr>
<td><img src="result-a-GT04.png"  width="320" ></td>
<td><img src="result-b-GT04.png"  width="320" ></td>
</tr>
<tr>
<td>estimated a (rescaled)</td>
<td>estimated b (rescaled)</td>
</tr>
</table>

The image is taken from the [alpha matting evaluation website](http://www.alphamatting.com) [3].

## References

```
[1] A. Levin, D. Lischinski, and Y. Weiss, “A closed-form solution to natural image matting,”
    Pattern Analysis and Machine Intelligence, IEEE Transactions on, vol. 30, no. 2,
    pp. 228–242, Feb 2008.
[2] Kaiming He, C. Rhemann, C. Rother, Xiaoou Tang, and Jian Sun, “A global sampling method
    for alpha matting,” in Computer Vision and Pattern Recognition (CVPR), 2011 IEEE Conference on,
	June 2011, pp. 2049–2056.
[3] Christoph Rhemann, Carsten Rother, Jue Wang, Margrit Gelautz, Pushmeet Kohli, and Pamela Rott,
    “Alpha matting evaluation website,” http://www.alphamatting.com.
```


