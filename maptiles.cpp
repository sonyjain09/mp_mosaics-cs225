/**
 * @file maptiles.cpp
 * Code for the maptiles function.
 */

#include <iostream>
#include <map>
#include "maptiles.h"
//#include "cs225/RGB_HSL.h"

using namespace std;


Point<3> convertToXYZ(LUVAPixel pixel) {
    return Point<3>( pixel.l, pixel.u, pixel.v );
}

MosaicCanvas* mapTiles(SourceImage const& theSource, vector<TileImage>& theTiles)
{
    int rows = theSource.getRows();
    int cols = theSource.getColumns();
    MosaicCanvas* output = new MosaicCanvas(rows, cols);
    std::vector<Point<3>> averages;
    map<Point<3>, TileImage*> m;
    for(TileImage& x : theTiles){
        Point<3> insert = convertToXYZ(x.getAverageColor());
        averages.push_back(insert);
        m[insert] = &x;
    }
    KDTree<3> tree(averages);
    for(int i = 0; i < rows; i++){
        for(int j = 0; j < cols; j++){
            Point<3> wanted = convertToXYZ(theSource.getRegionColor(i, j));
            Point<3> elem = tree.findNearestNeighbor(wanted);
            TileImage* neighbor = m[elem];
            output->setTile(i, j, neighbor);
        }
    }
    return output;
}