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

#include "eigen3/Eigen/SVD"
#include "iostream"
#include "opencv2/calib3d.hpp"
#include "opencv2/highgui/highgui.hpp"
#include "opencv2/stereo/stereo.hpp"
#include "opencv2/video/tracking.hpp"
#include "string"
#include <PointFeatureTracker.h>

using namespace GIFT;
using namespace Eigen;
using namespace std;
using namespace cv;

void PointFeatureTracker::processImage(const Mat& image) {
    this->trackFeatures(image);
    image.copyTo(this->previousImage);

    if (this->features.size() > this->featureSearchThreshold * this->maxFeatures)
        return;

    vector<Point2f> newPoints = this->detectNewFeatures(image);
    vector<Feature> newFeatures = this->createNewFeatures(image, newPoints);

    this->addNewFeatures(newFeatures);
}

vector<Feature> PointFeatureTracker::createNewFeatures(const Mat& image, const vector<Point2f>& newPoints) {
    vector<Feature> newFeatures;
    if (newPoints.empty())
        return newFeatures;

    for (int i = 0; i < newPoints.size(); ++i) {
        colorVec pointColor = {image.at<Vec3b>(newPoints[i]).val[0], image.at<Vec3b>(newPoints[i]).val[1],
            image.at<Vec3b>(newPoints[i]).val[2]};

        Feature lm(newPoints[i], cameraPtr, -1, pointColor);

        newFeatures.emplace_back(lm);
    }

    return newFeatures;
}

void PointFeatureTracker::trackFeatures(const Mat& image) {
    if (features.empty())
        return;

    vector<Point2f> oldPoints;
    for (const auto& feature : features) {
        oldPoints.emplace_back(feature.camCoordinates);
    }

    vector<Point2f> points;
    vector<uchar> status;
    vector<float> err;
    calcOpticalFlowPyrLK(previousImage, image, oldPoints, points, status, err, Size(winSize, winSize), maxLevel,
        TermCriteria((TermCriteria::COUNT) + (TermCriteria::EPS), 30, (0.01000000000000000021)));

    for (long int i = points.size() - 1; i >= 0; --i) {
        if (status[i] == 0 || err[i] >= maxError) {
            features.erase(features.begin() + i);
            continue;
        }

        if (!imageMask.empty()) {
            if (imageMask.at<uchar>(points[i]) == 0) {
                features.erase(features.begin() + i);
                continue;
            }
        }

        colorVec pointColor = {
            image.at<Vec3b>(points[i]).val[0], image.at<Vec3b>(points[i]).val[1], image.at<Vec3b>(points[i]).val[2]};
        features[i].update(points[i], pointColor);
    }
}

vector<Point2f> PointFeatureTracker::detectNewFeatures(const Mat& image) const {
    Mat imageGrey;
    cv::cvtColor(image, imageGrey, cv::COLOR_BGR2GRAY);

    vector<Point2f> proposedFeatures;
    goodFeaturesToTrack(imageGrey, proposedFeatures, maxFeatures, minHarrisQuality, featureDist, imageMask);
    vector<Point2f> newFeatures = this->removeDuplicateFeatures(proposedFeatures);

    return newFeatures;
}

vector<Point2f> PointFeatureTracker::removeDuplicateFeatures(const vector<Point2f>& proposedFeatures) const {
    vector<Point2f> newFeatures;
    for (const auto& proposedFeature : proposedFeatures) {
        bool useFlag = true;
        for (const auto& feature : this->features) {
            if (norm(proposedFeature - feature.camCoordinates) < featureDist) {
                useFlag = false;
                break;
            }
        }

        if (useFlag) {
            newFeatures.emplace_back(proposedFeature);
        }
    }
    return newFeatures;
}

void PointFeatureTracker::addNewFeatures(vector<Feature> newFeatures) {
    for (auto& lm : newFeatures) {
        if (features.size() >= maxFeatures)
            break;

        lm.idNumber = ++currentNumber;
        features.emplace_back(lm);
    }
}

void PointFeatureTracker::setMask(const Mat& mask, int cameraNumber) { imageMask = mask; }

Eigen::Matrix3T GIFT::skew_matrix(const Eigen::Vector3T& t) {
    Eigen::Matrix3T t_hat;
    t_hat << 0, -t(2), t(1), t(2), 0, -t(0), -t(1), t(0), 0;
    return t_hat;
}

Mat PointFeatureTracker::drawFeatureImage(const Scalar& color, const int pointSize, const int thickness) const {
    cv::Mat featureImage;
    this->previousImage.copyTo(featureImage);
    for (const auto& lm : this->features) {
        cv::circle(featureImage, lm.camCoordinates, pointSize, color, thickness);
    }
    return featureImage;
}

Mat PointFeatureTracker::drawFlowImage(
    const Scalar& featureColor, const Scalar& flowColor, const int pointSize, const int thickness) const {
    Mat flowImage = drawFeatureImage(featureColor, pointSize, thickness);
    for (const auto& lm : this->features) {
        Point2f p1 = lm.camCoordinates;
        Point2f p0 = p1 - Point2f(lm.opticalFlowRaw.x(), lm.opticalFlowRaw.y());
        line(flowImage, p0, p1, flowColor, thickness);
    }
    return flowImage;
}

Mat PointFeatureTracker::drawFlow(
    const Scalar& featureColor, const Scalar& flowColor, const int pointSize, const int thickness) const {
    Mat flow(this->previousImage.size(), CV_8UC3);
    flow.setTo(0);

    for (const auto& lm : this->features) {
        Point2f p1 = lm.camCoordinates;
        Point2f p0 = p1 - Point2f(lm.opticalFlowRaw.x(), lm.opticalFlowRaw.y());
        circle(flow, p1, pointSize, featureColor, thickness);
        line(flow, p0, p1, flowColor, thickness);
    }
    return flow;
}

/*
void PointFeatureTracker::computeLandmarkPositions() {
    if (mode == TrackerMode::MONO) return;
    for (auto & lm : features) {
        if (mode == TrackerMode::STEREO) {
            lm.position = this->solveStereo(lm.camCoordinatesNorm()[0], lm.camCoordinatesNorm()[1]);
        }
    }
}

Vector3T PointFeatureTracker::solveStereo(const Point2f& leftKp, const Point2f& rightKp) const {
    assert(leftKp.x > rightKp.x);
    Vector3T position;
    position << leftKp.x, leftKp.y, 1;

    ftype scale = stereoBaseline / (leftKp.x - rightKp.x);
    position = scale*position;

    return position;
}

Vector3T PointFeatureTracker::solveMultiView(const vector<Point2f> imageCoordinatesNorm) const {
    int camNum = cameras.size();
    assert(imageCoordinatesNorm.size() == cameras.size());
    assert(camNum >= 2);

    MatrixXT solutionMat(3*camNum, 4);

    for (int i=0; i<camNum; ++i) {
        Vector3T imageCoord;
        imageCoord << imageCoordinatesNorm[i].x, imageCoordinatesNorm[i].y, 1;
        solutionMat.block<3,4>(3*i,0) = skew_matrix(imageCoord) * cameras[i].P;
    }

    JacobiSVD<MatrixXT> svd(solutionMat, ComputeThinU | ComputeFullV);
    Vector4T positionHomogeneous = svd.matrixV().block<4,1>(0,3);

    Vector3T position = positionHomogeneous.block<3,1>(0,0) / positionHomogeneous(3);
    return position;
}

vector<vector<Point2f>> PointFeatureTracker::detectNewStereoFeatures(const cv::Mat &imageLeft, const cv::Mat
&imageRight) const { vector<vector<Point2f>> newFeatures(2);

    // Obtain left features
    Mat imageGreyLeft;
    cv::cvtColor(imageLeft, imageGreyLeft, cv::COLOR_BGR2GRAY);
    goodFeaturesToTrack(imageGreyLeft, newFeatures[0], 2*maxFeatures, 0.001, featureDist, imageMasks[0]);
    newFeatures[0] = this->removeDuplicateFeatures(newFeatures[0]);

    // Track features to the right image
    vector<uchar> status;
    vector<float> err;
    calcOpticalFlowPyrLK(imageLeft, imageRight, newFeatures[0], newFeatures[1], status, err);

    assert(newFeatures[0].size() == newFeatures[1].size());

    for (long int i=newFeatures[0].size()-1; i >= 0; --i) {
        bool eraseCondition = (status[i] == 0);
        if (!imageMasks[0].empty()) eraseCondition |= imageMasks[0].at<uchar>(newFeatures[0][i]);
        if (!imageMasks[1].empty()) eraseCondition |= imageMasks[1].at<uchar>(newFeatures[1][i]);

        if (eraseCondition) {
            newFeatures[0].erase(newFeatures[0].begin() +i);
            newFeatures[1].erase(newFeatures[1].begin() +i);
        }
    }

    return newFeatures;

}
*/
