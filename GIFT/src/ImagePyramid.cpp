/* 
    This file is part of GIFT.

    GIFT is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    GIFT is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with GIFT.  If not, see <http://www.gnu.org/licenses/>.
*/

#include "ImagePyramid.h"
#include "opencv2/imgproc/imgproc.hpp"
#include <algorithm>

ImagePyramid::ImagePyramid(const cv::Mat& image, const int& numLevels) {
    assert(numLevels > 0);
    levels.resize(numLevels);
    levels[0] = image;
    for (int i=1; i<numLevels; ++i) {
        cv::pyrDown(levels[i-1], levels[i]);
    }
}

ImageWithGradient::ImageWithGradient(const cv::Mat& image) {
    this->image = image;
    cv::Sobel(image, this->gradientX, CV_32F, 1, 0, -1, 1.0, 0.0, cv::BORDER_REPLICATE);
    cv::Sobel(image, this->gradientY, CV_32F, 0, 1, -1, 1.0, 0.0, cv::BORDER_REPLICATE);
}

ImageWithGradientPyramid::ImageWithGradientPyramid(const cv::Mat& image, const int& numLevels) {
    assert(numLevels > 0);
    levels.resize(numLevels);
    levels[0] = ImageWithGradient(image);
    for (int i=1; i<numLevels; ++i) {
        cv::Mat temp;
        cv::pyrDown(levels[i-1].image, temp);
        levels[i] = ImageWithGradient(temp);
    }
}

PyramidPatch extractPyramidPatch(const cv::Point2f& point, const cv::Size& sze, const ImageWithGradientPyramid& pyr) {
    int numLevels = pyr.levels.size();
    PyramidPatch patch;
    patch.baseCentre = Vector2T(point.x, point.y);
    patch.vecImage.resize(numLevels);
    patch.vecDifferential.resize(numLevels);
    patch.rows = sze.height; patch.cols = sze.width;
    for (int lv=0; lv<numLevels; ++lv) {
        Mat tempI, tempX, tempY;
        getRectSubPix(pyr.levels[lv].image,     sze, point*pow(2,-lv), tempI);
        getRectSubPix(pyr.levels[lv].gradientX, sze, point*pow(2,-lv), tempX);
        getRectSubPix(pyr.levels[lv].gradientY, sze, point*pow(2,-lv), tempY);
        
        patch.vecImage[lv] = vectoriseImage(tempI);
        Matrix<ftype, Dynamic, 2> temp(sze.area(),2);
        temp.block(0,0,sze.area(),1) = vectoriseImage(tempX);
        temp.block(0,1,sze.area(),1) = vectoriseImage(tempY);
        patch.vecDifferential[lv] = temp;
    }
    return patch;
}

vector<PyramidPatch> extractPyramidPatches(const vector<cv::Point2f>& points, const cv::Mat& image, const cv::Size& sze, const int& numLevels) {
    ImageWithGradientPyramid pyr(image, numLevels);
    auto patchLambda = [pyr, sze](const cv::Point2f& point) { return extractPyramidPatch(point, sze, pyr); };
    vector<PyramidPatch> patches;
    transform(points.begin(), points.end(), patches.begin(), patchLambda);
    return patches;
}

ImagePatch getPatchAtLevel(const PyramidPatch& pyrPatch, const int lv) {
    ImagePatch patch;
    patch.vecImage = pyrPatch.vecImage[lv];
    patch.vecDifferential = pyrPatch.vecDifferential[lv];
    patch.centre = pyrPatch.baseCentre * pow(2,-lv);
    patch.rows = pyrPatch.rows; patch.cols = pyrPatch.cols;
    return patch;
}

VectorXT vectoriseImage(const Mat& image) {
    // We work row by row
    const int rows = image.rows;
    const int cols = image.cols;
    VectorXT vecImage(rows*cols);
    for (int y=0; y<rows; ++y) {
    for (int x=0; x<cols; ++x) {
        switch ( image.depth() ) {
            case CV_8U : vecImage(x+y*rows) = (ftype) image.at<uchar>(Point2i(x,y)); break;
            case CV_16S : vecImage(x+y*rows) = (ftype) image.at<short>(Point2i(x,y)); break;
            case CV_32F : vecImage(x+y*rows) = (ftype) image.at<float>(Point2i(x,y)); break;
            case CV_64F : vecImage(x+y*rows) = (ftype) image.at<double>(Point2i(x,y)); break;
        }
    }
    }
    return vecImage;
}