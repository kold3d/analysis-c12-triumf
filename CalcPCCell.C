typedef struct Vect3d {
  Float_t x;
  Float_t y;
  Float_t z;
} Point3d;

class Vec3d {
public:
  Vec3d(Point3d start, Point3d stop) {
    mag = sqrt((start.x-stop.x)*(start.x-stop.x)+
	       (start.y-stop.y)*(start.y-stop.y)+
	       (start.z-stop.z)*(start.z-stop.z));
    x=(start.x-stop.x)/mag;
    y=(start.y-stop.y)/mag;
    z=(start.z-stop.z)/mag;
  };
  Point3d PointIntersectXPlane(Float_t xp) {
   Point3d point = {xp,xp/x*y,xp/x*z};
   return point;
  };
  Point3d PointIntersectYPlane(Float_t yp) {
    Point3d point = {yp/y*x,yp,yp/y*z};
    return point;
  };
  Point3d PointIntersectZPlane(Float_t zp) {
    Point3d point = {zp/z*x,zp/z*y,zp};
    return point;
  };
  Float_t mag;
  Float_t x;
  Float_t y;
  Float_t z;
};

typedef std::pair<Float_t,Float_t> Boundary;
		  
std::pair<Int_t,Float_t> CalcPCCell(Float_t x1, Float_t y1, Float_t depth) {
		    
  Float_t SiToBack  = 28.5-5.08;
  Float_t SiToFront = 28.5+5.08;
  Float_t z2 = 513.-depth;
  
  
  Boundary cell[5];
  cell[0] = Boundary( -25.4,-15.24);
  cell[1] = Boundary(-15.24,-5.08);
  cell[2] = Boundary( -5.08,  5.08);
  cell[3] = Boundary(  5.08, 15.24);
  cell[4] = Boundary( 15.24, 25.4 );
  
  Point3d siPoint = {x1,y1,0.};
  Point3d chamberPoint = {0.,0.,z2};
  Vec3d vec(siPoint,chamberPoint);
  
  Char_t entranceCell = -1;
  Point3d entrancePoint = vec.PointIntersectZPlane(SiToBack);
  for(UChar_t i = 0;i<5;i++) {
    if(cell[i].first < entrancePoint.y && entrancePoint.y <= cell[i].second) {
      entranceCell = i;
      break;
    }
  }

  Char_t exitCell = -1;
  Point3d exitPoint = vec.PointIntersectZPlane(SiToFront);
  for(UChar_t i = 0;i<5;i++) {
    if(cell[i].first < exitPoint.y && exitPoint.y <= cell[i].second) {
      exitCell = i;
      break;
    }
  }

  std::pair<Int_t,Float_t> returnValues;
  if(entranceCell != exitCell) {
    
  } else {
    returnValues.first = entranceCell;
    
  }
  
  return returnValues;
}
