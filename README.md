# RBF-based_correlated_random_field_generator

![Preview](https://numericalenvironmental.files.wordpress.com/2019/06/3d-field-image2.jpg?w=1632)

This is a streamlined python 3 script for generating spatially-correlated random fields in 2-D or 3-D using a radial basis function interpolator. It furnishes what I think are significant improvements in simplicity as well as function compared to an earlier effort to do the same (link to 1st posting). This script requires pandas as well as several scipy tools. See my blog post (https://numericalenvironmental.wordpress.com/2019/06/01/an-improved-3-d-correlated-random-field-generator-in-python/) for a brief discussion.

The following text input files are required:

* seeds.csv - seed points and values; columns must be filled in the all four fields [x, y, z, v], even for 2-D problems
* params.txt - miscellaneous model setup parameters (e.g., gridding, anisotropy, search radius, etc.)

I'd appreciate hearing back from you if you find the code useful. Questions or comments are welcome at walt.mcnab@gmail.com.

THIS CODE/SOFTWARE IS PROVIDED IN SOURCE OR BINARY FORM "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

