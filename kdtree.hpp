/**
 * @file kdtree.cpp
 * Implementation of KDTree class.
 */

#include <utility>
#include <algorithm>
#include <iostream>
#include <cmath>
#include <cstdlib>
#include <math.h>

using namespace std;

template <int Dim>
bool KDTree<Dim>::smallerDimVal(const Point<Dim>& first, const Point<Dim>& second, int curDim) const{
  if (first[curDim] == second[curDim]) return first < second;
  return first[curDim] < second[curDim];
}

template <int Dim>
double KDTree<Dim>::getDist(const Point<Dim>& check, const Point<Dim>& target) const{
  double dist = 0;
  for(int i =  0; i < Dim; i++) dist += pow((target[i]-check[i]), 2);
  return dist;
}

template <int Dim>
bool KDTree<Dim>::shouldReplace(const Point<Dim>& target, const Point<Dim>& currentBest, const Point<Dim>& potential) const{
  double pDist = getDist(potential, target);
  double cDist = getDist(currentBest, target);
  if(pDist == cDist) return potential < currentBest;
  return pDist < cDist;
}

template <int Dim>
void KDTree<Dim>::swap(vector<Point<Dim>>& newPoints, int one, int two){
  Point<Dim> temp = newPoints[one];
  newPoints[one] = newPoints[two];
  newPoints[two] = temp;
}

template <int Dim>
int KDTree<Dim>::partition(vector<Point<Dim>>& newPoints, int left, int right, int pivot, int curDim){
  Point<Dim> val = newPoints[pivot];
  swap(newPoints, pivot, right);
  int iter = left;
  for(int i = left; i < right; i++){
    if(smallerDimVal(newPoints[i],val, curDim)){
      swap(newPoints, iter, i);
      iter++;
    }
  }
  swap(newPoints, iter, right);
  return iter;
}

template <int Dim>
Point<Dim>& KDTree<Dim>::select(vector<Point<Dim>>& newPoints, int left, int right, int target, int curDim){
  if(left == right) return newPoints[left];
  int pivot = partition(newPoints, left, right, target, curDim%Dim);
  if(pivot == target) return newPoints[target];
  else if(pivot > target) return select(newPoints, left, pivot-1, target, curDim%Dim);
  else return select(newPoints, pivot+1, right, target, curDim%Dim);
}

template <int Dim>
typename KDTree<Dim>::KDTreeNode* KDTree<Dim>::treeBuilder(vector<Point<Dim>>& newPoints, int left, int right, int curDim){
  if(left > right) return nullptr;
  int medIdx = (left+right)/2;
  Point<Dim> median = select(newPoints, left, right, medIdx, curDim%Dim);
  KDTreeNode* subroot = new KDTreeNode(median);
  size++;
  if (left != right) {
    subroot->left = treeBuilder(newPoints, left, medIdx - 1, (curDim + 1)%Dim);
    subroot->right = treeBuilder(newPoints, medIdx + 1, right, (curDim + 1)%Dim);
  }
  return subroot;
} 

template <int Dim>
KDTree<Dim>::KDTree(const vector<Point<Dim>>& newPoints){
  size = 0;
  if (newPoints.size() == 0) {
    root = nullptr;
    return;
  }
  vector<Point<Dim>> copy = newPoints;
  root = treeBuilder(copy, 0, newPoints.size()-1, 0);
}

template <int Dim>
KDTree<Dim>::KDTree(const KDTree<Dim>& other) {
  clear(root);
  copy(this->root, other.root);
  size = other.size;
}

template <int Dim>
const KDTree<Dim>& KDTree<Dim>::operator=(const KDTree<Dim>& rhs) {
  clear(root);
  copy(this->root, rhs.root);
  size = rhs.size;
  return *this;
}

template <int Dim>
KDTree<Dim>::~KDTree() {
  clear(root);
}

template <int Dim>
void KDTree<Dim>::clear(KDTreeNode* subroot) {
  if(subroot == NULL) return;
  clear(subroot->left);
  clear(subroot->right);
  delete subroot;
}

template <int Dim>
void KDTree<Dim>::copy(KDTreeNode* self, KDTreeNode* other){
  if(other == NULL) return;
  self->point = other->point;
  copy(self->left, other->left);
	copy(self->right, other->right);
}

template <int Dim>
Point<Dim> KDTree<Dim>::findNearestNeighbor(const Point<Dim>& query) const{
  return neighborHelper(root, query, 0);
}

template <int Dim>
Point<Dim> KDTree<Dim>::neighborHelper(KDTreeNode* subroot, const Point<Dim>& query, int curDim) const {
  bool left = false;
	Point<Dim> nearest = subroot->point;
	if (subroot->left == NULL and subroot->right == NULL) return nearest;
	if (smallerDimVal(query, nearest, curDim)) {
    if(subroot->left != NULL) nearest = neighborHelper(subroot->left, query, (curDim+1) % Dim);
    left = true;
	}
	else if(subroot->right != NULL) nearest = neighborHelper(subroot->right, query, (curDim + 1) % Dim);
  else return nearest;
	if (shouldReplace(query, nearest, subroot->point)) nearest = subroot->point;
  double radius = getDist(query, nearest);
  double splitDist = pow((subroot->point[curDim]-query[curDim]),2);
  if (radius >= splitDist) {
    Point<Dim> other;
    if(left and subroot->right != NULL) other = neighborHelper(subroot->right, query, (curDim + 1) % Dim);
    else if(subroot->left != NULL) other = neighborHelper(subroot->left, query, (curDim + 1) % Dim);
    else return nearest;
    if(shouldReplace(query, nearest, other)) nearest = other;
  }
	return nearest;
} 