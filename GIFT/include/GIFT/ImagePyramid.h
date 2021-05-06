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

#pragma once

#include "eigen3/Eigen/Core"
#include "ftype.h"
#include "opencv2/core/core.hpp"
#include <vector>

namespace GIFT {

struct ImagePyramid {
    std::vector<cv::Mat> levels;
    ImagePyramid(){};
    ImagePyramid(const cv::Mat& image, const int& numLevels);
};

struct ImageWithGradient {
    cv::Mat image;
    cv::Mat gradientX;
    cv::Mat gradientY;
    ImageWithGradient(){};
    ImageWithGradient(const cv::Mat& image);
};

struct ImageWithGradientPyramid {
    std::vector<ImageWithGradient> levels;
    ImageWithGradientPyramid(){};
    ImageWithGradientPyramid(const cv::Mat& image, const int& numLevels);
    ImageWithGradientPyramid(const ImagePyramid& imagePyr);
};

struct ImagePatch {
    ImageWithGradient imageWithGrad;
    Eigen::Vector2T centre;
    int rows() const { return imageWithGrad.image.rows; };
    int cols() const { return imageWithGrad.image.cols; };
    int area() const { return rows() * cols(); }
    ftype at(int row, int col) const;
    Eigen::Matrix<ftype, 1, 2> differential(int row, int col) const;

    Eigen::VectorXT imageVector() const;
    Eigen::Matrix<ftype, Eigen::Dynamic, 2> imageVectorDifferential() const;
};

struct PyramidPatch {
    std::vector<ImagePatch> levels;

    int rows(const int& lv = 0) const { return levels[lv].rows(); };
    int cols(const int& lv = 0) const { return levels[lv].cols(); };
    ftype at(int row, int col, int lv = 0) const;
    int totalPixelCount() const;
    Eigen::Vector2T centre(const int& lv = 0) const { return levels[lv].centre; };
    Eigen::VectorXT pyramidVector() const;
    Eigen::Matrix<ftype, Eigen::Dynamic, 2> pyramidVectorDifferential() const;
};

PyramidPatch extractPyramidPatch(const cv::Point2f& point, const cv::Size& sze, const ImageWithGradientPyramid& pyr);
ImagePatch extractImagePatch(const cv::Point2f& point, const cv::Size& sze, const ImageWithGradient& imageWithGrad);
PyramidPatch extractPyramidPatch(
    const cv::Point2f& point, const std::vector<cv::Size>& sizes, const ImageWithGradientPyramid& pyr);
std::vector<PyramidPatch> extractPyramidPatches(
    const std::vector<cv::Point2f>& points, const cv::Mat& image, const cv::Size& sze, const int& numLevels);
Eigen::VectorXT vectoriseImage(const cv::Mat& image);
ftype pixelValue(const cv::Mat& image, const int& row, const int& col);

} // namespace GIFT