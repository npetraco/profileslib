#include <Rcpp.h>
#include <iostream>
#include <fstream>
#include <sys/stat.h>
#include <vector>

using namespace std;
using namespace Rcpp;

//Helper function prototypes:
int * index2sub(int, int);


// [[Rcpp::export]]
List readSurFile(std::string path) {

  struct stat results; //General file info
  const char* path_ptr;
  path_ptr=path.c_str();

  //Gather some general file info:
  stat(path_ptr, &results);

  try{
    //Open .sur file
    ifstream surFileStream(path_ptr, ios::in | ios::binary);

    int streamIdx = 0;

    // 1
    //Read in info from main header:
    char bcode[12]="";
    surFileStream.seekg(streamIdx);
    surFileStream.read((char *)&bcode, sizeof(bcode));
    streamIdx+=sizeof(bcode);
    //cout<<bcode<<endl;

    // 2
    unsigned short format;
    surFileStream.seekg(streamIdx);
    surFileStream.read((char *)&format, sizeof(format));
    streamIdx+=sizeof(format);
    //cout<<format<<endl;

    // 3
    unsigned short numObjects;
    surFileStream.seekg(streamIdx);
    surFileStream.read((char *)&numObjects, sizeof(numObjects));
    streamIdx+=sizeof(numObjects);
    //cout<<numObjects<<endl;

    // 4
    unsigned short verNum; // Should be 1
    surFileStream.seekg(streamIdx);
    surFileStream.read((char *)&verNum, sizeof(verNum));
    streamIdx+=sizeof(verNum);
    //cout<<verNum<<endl;

    // 5
    unsigned short studiableType; // 1 for a profile
                                  // 2 for a surface
                                  // 3 for binary image
                                  // 4 for series of profiles
									                // 5 for series of surfaces
                                  // 6 for meridian discs
    surFileStream.seekg(streamIdx);
    surFileStream.read((char *)&studiableType, sizeof(studiableType));
    streamIdx+=sizeof(studiableType);
    //cout<<studiableType<<endl;

    // 6
	  char nameObject[30]="";
    surFileStream.seekg(streamIdx);
    surFileStream.read((char *)&nameObject, sizeof(nameObject));
    streamIdx+=sizeof(nameObject);
    //cout<<nameObject<<endl;

    // 7
    char nameOperator[30]="";
    surFileStream.seekg(streamIdx);
    surFileStream.read((char *)&nameOperator, sizeof(nameOperator));
    streamIdx+=sizeof(nameOperator);
    //cout<<nameOperator<<endl;

    // 8
    unsigned short unused1;
    surFileStream.seekg(streamIdx);
    surFileStream.read((char *)&unused1, sizeof(unused1));
    streamIdx+=sizeof(unused1);
    //cout<<unused1<<endl;

    // 9
    unsigned short unused2;
    surFileStream.seekg(streamIdx);
    surFileStream.read((char *)&unused2, sizeof(unused2));
    streamIdx+=sizeof(unused2);
    //cout<<unused2<<endl;

    // 10
    unsigned short unused3;
    surFileStream.seekg(streamIdx);
    surFileStream.read((char *)&unused3, sizeof(unused3));
    streamIdx+=sizeof(unused3);
    //cout<<unused3<<endl;

    // 11
    unsigned short non_meas_pts_flg;
    surFileStream.seekg(streamIdx);
    surFileStream.read((char *)&non_meas_pts_flg, sizeof(non_meas_pts_flg));
    streamIdx+=sizeof(non_meas_pts_flg);
    //cout<<non_meas_pts_flg<<endl;

    // 12
    unsigned short absZAxisFlag;
    surFileStream.seekg(streamIdx);
    surFileStream.read((char *)&absZAxisFlag, sizeof(absZAxisFlag));
    streamIdx+=sizeof(absZAxisFlag);
    //cout<<absZAxisFlag<<endl;

    // 13
    char reserved[8]=""; //Should be a byte, but Rcpp seems to have trouble with a template of unsigned chars
    surFileStream.seekg(streamIdx);
    surFileStream.read((char *)&reserved, sizeof(reserved));
    streamIdx+=sizeof(reserved);
    //cout<<reserved<<endl;

    // 14
    unsigned short num_bits_pt; // Number of bits per point
    surFileStream.seekg(streamIdx);
    surFileStream.read((char *)&num_bits_pt, sizeof(num_bits_pt));
    streamIdx+=sizeof(num_bits_pt);
    //cout<<num_bits_pt<<endl;

    // 15
	  int minPt;
    surFileStream.seekg(streamIdx);
    surFileStream.read((char *)&minPt, sizeof(minPt));
    streamIdx+=sizeof(minPt);
    //cout<<minPt<<endl;

    // 16
    int maxPt;
    surFileStream.seekg(streamIdx);
    surFileStream.read((char *)&maxPt, sizeof(maxPt));
    streamIdx+=sizeof(maxPt);
    //cout<<maxPt<<endl;

    // 17
    int num_pts_profile; // Width of surface
    surFileStream.seekg(streamIdx);
    surFileStream.read((char *)&num_pts_profile, sizeof(num_pts_profile));
    streamIdx+=sizeof(num_pts_profile);
    //cout<<num_pts_profile<<endl;

    // 18
  	int num_profiles; // Length of surface
    surFileStream.seekg(streamIdx);
    surFileStream.read((char *)&num_profiles, sizeof(num_profiles));
    streamIdx+=sizeof(num_profiles);
    //cout<<num_profiles<<endl;

    // 19
	  int totalNumPts;
    surFileStream.seekg(streamIdx);
    surFileStream.read((char *)&totalNumPts, sizeof(totalNumPts));
    streamIdx+=sizeof(totalNumPts);
    //cout<<totalNumPts<<endl;

    // 20
    float x_inc; // Sampling interval in x-direction.
    surFileStream.seekg(streamIdx);
    surFileStream.read((char *)&x_inc, sizeof(x_inc));
    streamIdx+=sizeof(x_inc);
    //cout<<x_inc<<endl;

    // 21
	  float y_inc; // Sampling interval in y-direction.
    surFileStream.seekg(streamIdx);
    surFileStream.read((char *)&y_inc, sizeof(y_inc));
    streamIdx+=sizeof(y_inc);
    //cout<<y_inc<<endl;

    // 22
    float z_inc; // Sampling interval in y-direction.
    surFileStream.seekg(streamIdx);
    surFileStream.read((char *)&z_inc, sizeof(z_inc));
    streamIdx+=sizeof(z_inc);
    //cout<<z_inc<<endl;

    // 23
    char xname[16]="";
    surFileStream.seekg(streamIdx);
    surFileStream.read((char *)&xname, sizeof(xname));
    streamIdx+=sizeof(xname);
    //cout<<xname<<endl;

    // 24
    char yname[16]="";
    surFileStream.seekg(streamIdx);
    surFileStream.read((char *)&yname, sizeof(yname));
    streamIdx+=sizeof(yname);
    //cout<<yname<<endl;

    // 25
    char zname[16]="";
    surFileStream.seekg(streamIdx);
    surFileStream.read((char *)&zname, sizeof(zname));
    streamIdx+=sizeof(zname);
    //cout<<zname<<endl;

    // 26
    char xstepunits[16]="";
    surFileStream.seekg(streamIdx);
    surFileStream.read((char *)&xstepunits, sizeof(xstepunits));
    streamIdx+=sizeof(xstepunits);
    //cout<<xstepunits<<endl;

    // 27
    char ystepunits[16]="";
    surFileStream.seekg(streamIdx);
    surFileStream.read((char *)&ystepunits, sizeof(ystepunits));
    streamIdx+=sizeof(ystepunits);
    //cout<<ystepunits<<endl;

    // 28
    char zstepunits[16]="";
    surFileStream.seekg(streamIdx);
    surFileStream.read((char *)&zstepunits, sizeof(zstepunits));
    streamIdx+=sizeof(zstepunits);
    //cout<<zstepunits<<endl;

    // 29
    char xunits[16]="";
    surFileStream.seekg(streamIdx);
    surFileStream.read((char *)&xunits, sizeof(xunits));
    streamIdx+=sizeof(xunits);
    //cout<<xunits<<endl;

    // 30
    char yunits[16]="";
    surFileStream.seekg(streamIdx);
    surFileStream.read((char *)&yunits, sizeof(yunits));
    streamIdx+=sizeof(yunits);
    //cout<<yunits<<endl;

    // 31
    char zunits[16]="";
    surFileStream.seekg(streamIdx);
    surFileStream.read((char *)&zunits, sizeof(zunits));
    streamIdx+=sizeof(zunits);
    //cout<<zunits<<endl;

    // 32
    float xratio;
    surFileStream.seekg(streamIdx);
    surFileStream.read((char *)&xratio, sizeof(xratio));
    streamIdx+=sizeof(xratio);
    //cout<<xratio<<endl;

    // 33
    float yratio;
    surFileStream.seekg(streamIdx);
    surFileStream.read((char *)&yratio, sizeof(yratio));
    streamIdx+=sizeof(yratio);
    //cout<<yratio<<endl;

    // 34
    float zratio;
    surFileStream.seekg(streamIdx);
    surFileStream.read((char *)&zratio, sizeof(zratio));
    streamIdx+=sizeof(zratio);
    //cout<<zratio<<endl;

    // 35
    short replica;
    surFileStream.seekg(streamIdx);
    surFileStream.read((char *)&replica, sizeof(replica));
    streamIdx+=sizeof(replica);
    //cout<<replica<<endl;

    // 36
	  short inverted;
    surFileStream.seekg(streamIdx);
    surFileStream.read((char *)&inverted, sizeof(inverted));
    streamIdx+=sizeof(inverted);
    //cout<<inverted<<endl;

    // 37
	  short leveled;
    surFileStream.seekg(streamIdx);
    surFileStream.read((char *)&leveled, sizeof(leveled));
    streamIdx+=sizeof(leveled);
    //cout<<leveled<<endl;

    // 38
	  float unused4;
    surFileStream.seekg(streamIdx);
    surFileStream.read((char *)&unused4, sizeof(unused4));
    streamIdx+=sizeof(unused4);
    //cout<<unused4<<endl;

    // 39
    float unused5;
    surFileStream.seekg(streamIdx);
    surFileStream.read((char *)&unused5, sizeof(unused5));
    streamIdx+=sizeof(unused5);
    //cout<<unused5<<endl;

    // 40
    float unused6;
    surFileStream.seekg(streamIdx);
    surFileStream.read((char *)&unused6, sizeof(unused6));
    streamIdx+=sizeof(unused6);
    //cout<<unused6<<endl;

    // 41
	  short secs;
    surFileStream.seekg(streamIdx);
    surFileStream.read((char *)&secs, sizeof(secs));
    streamIdx+=sizeof(secs);
    //cout<<secs<<endl;

    // 42
	  short mins;
    surFileStream.seekg(streamIdx);
    surFileStream.read((char *)&mins, sizeof(mins));
    streamIdx+=sizeof(mins);
    //cout<<mins<<endl;

    // 43
	  short hrs;
    surFileStream.seekg(streamIdx);
    surFileStream.read((char *)&hrs, sizeof(hrs));
    streamIdx+=sizeof(hrs);
    //cout<<hrs<<endl;

    // 44
	  short day;
    surFileStream.seekg(streamIdx);
    surFileStream.read((char *)&day, sizeof(day));
    streamIdx+=sizeof(day);
    //cout<<day<<endl;

    // 45
	  short month;
    surFileStream.seekg(streamIdx);
    surFileStream.read((char *)&month, sizeof(month));
    streamIdx+=sizeof(month);
    //cout<<month<<endl;

    // 46
	  short year;
    surFileStream.seekg(streamIdx);
    surFileStream.read((char *)&year, sizeof(year));
    streamIdx+=sizeof(year);
    //cout<<year<<endl;

    // 47
	  short wkday;
    surFileStream.seekg(streamIdx);
    surFileStream.read((char *)&wkday, sizeof(wkday));
    streamIdx+=sizeof(wkday);
    //cout<<wkday<<endl;

    // 48
    //cout<<streamIdx<<endl;
	  float duration;
    surFileStream.seekg(streamIdx);
    surFileStream.read((char *)&duration, sizeof(duration));
    streamIdx+=sizeof(duration);
    //cout<<duration<<endl;

    // 49
	  short unused7;
    surFileStream.seekg(streamIdx);
    surFileStream.read((char *)&unused7, sizeof(unused7));
    streamIdx+=sizeof(unused7);
    //cout<<unused7<<endl;

    //  50
    short unused8;
    surFileStream.seekg(streamIdx);
    surFileStream.read((char *)&unused8, sizeof(unused8));
    streamIdx+=sizeof(unused8);
    //cout<<unused8<<endl;

    // 51
    float unused9;
    surFileStream.seekg(streamIdx);
    surFileStream.read((char *)&unused9, sizeof(unused9));
    streamIdx+=sizeof(unused9);
    //cout<<unused9<<endl;

    //52
    short unused10;
    surFileStream.seekg(streamIdx);
    surFileStream.read((char *)&unused10, sizeof(unused10));
    streamIdx+=sizeof(unused10);
    //cout<<unused10<<endl;

    // 53
	  short leng_comment_zone; //
    surFileStream.seekg(streamIdx);
    surFileStream.read((char *)&leng_comment_zone, sizeof(leng_comment_zone));
    streamIdx+=sizeof(leng_comment_zone);
    //cout<<leng_comment_zone<<endl;

    // 54
    //cout<<streamIdx<<endl;
	  short leng_private_zone; //
    surFileStream.seekg(streamIdx);
    surFileStream.read((char *)&leng_private_zone, sizeof(leng_private_zone));
    streamIdx+=sizeof(leng_private_zone);
    //cout<<leng_private_zone<<endl;

    //55
	  char freeZone[128];   //Should be a byte, but Rcpp seems to have trouble with a template of unsigned chars
    surFileStream.seekg(streamIdx);
    surFileStream.read((char *)&freeZone, sizeof(freeZone));
    streamIdx+=sizeof(freeZone);
    //cout<<freeZone<<endl;

    // 56
	  float xoff;
    surFileStream.seekg(streamIdx);
    surFileStream.read((char *)&xoff, sizeof(xoff));
    streamIdx+=sizeof(xoff);
    //cout<<xoff<<endl;

    // 57
    float yoff;
    surFileStream.seekg(streamIdx);
    surFileStream.read((char *)&yoff, sizeof(yoff));
    streamIdx+=sizeof(yoff);
    //cout<<yoff<<endl;

    // 58
    float zoff;
    surFileStream.seekg(streamIdx);
    surFileStream.read((char *)&zoff, sizeof(zoff));
    streamIdx+=sizeof(zoff);
    //cout<<zoff<<endl;

    // 59
    float Tspacing;
    surFileStream.seekg(streamIdx);
    surFileStream.read((char *)&Tspacing, sizeof(Tspacing));
    streamIdx+=sizeof(Tspacing);
    //cout<<Tspacing<<endl;

    // 60
	  float Toff;
    surFileStream.seekg(streamIdx);
    surFileStream.read((char *)&Toff, sizeof(Toff));
    streamIdx+=sizeof(Toff);
    //cout<<Toff<<endl;

    // 61
	  char TName[13];
    surFileStream.seekg(streamIdx);
    surFileStream.read((char *)&TName, sizeof(TName));
    streamIdx+=sizeof(TName);
    //cout<<TName<<endl;

    // 62
	  char TUnit[13];
    surFileStream.seekg(streamIdx);
    surFileStream.read((char *)&TUnit, sizeof(TUnit));
    streamIdx+=sizeof(TUnit);
    //cout<<TUnit<<endl;

    //Read in the byte array of surface points:
    int offsetToImage = (512+leng_comment_zone+leng_private_zone); //Handy to have for external programs

    //Read in the surface points:
    surFileStream.seekg(offsetToImage);
    //cout<<surFileStream.tellg()<<endl;
    NumericVector real_surface(totalNumPts);
    //NumericMatrix real_surface(num_profiles,num_pts_profile);

    if(num_bits_pt==16) {
      std::vector<short> vs;
      vs.resize( totalNumPts );
      surFileStream.read( (char*)( &vs[0]), totalNumPts * sizeof(short) );
      for(int i=0; i<totalNumPts; i++) {
        real_surface(i) =  1000 * z_inc * vs.at(i);
        //int * subscr = index2sub(i, num_pts_profile);
        //real_surface(subscr[0],subscr[1]) =  1000 * z_inc * vs.at(i);
        //int r = ((i-1) % num_pts_profile) + 1;
        //int c = floor((i-1) / num_pts_profile) + 1;
        //real_surface(r,c) =  1000 * z_inc * vs.at(i);

      }
    }

    if(num_bits_pt==32) {
      std::vector<int> vs;
      vs.resize( totalNumPts );
      surFileStream.read( (char*)( &vs[0]), totalNumPts * sizeof(int) );
      for(int i=0; i<totalNumPts; i++) {
        real_surface(i) =  1000 * z_inc * vs.at(i);
        //int * subscr = index2sub(i, num_pts_profile);
        //real_surface(subscr[0],subscr[1]) =  1000 * z_inc * vs.at(i);
      }
    }

    surFileStream.close();

    //Have to break this up into chuncks of <20. There is a limit on List length in the current Rcpp.
    //Is there a better way to do this with wrap????
    List info1;
    info1 = List::create(
         Rcpp::Named("code") = bcode,
         Rcpp::Named("format") = format,
         Rcpp::Named("num.obj")= numObjects,
         Rcpp::Named("ver.num")= verNum,
         Rcpp::Named("stud.typ")= studiableType,
         Rcpp::Named("name.obj")= nameObject,
         Rcpp::Named("name.op")= nameOperator,
         Rcpp::Named("unused.int1")= unused1,
         Rcpp::Named("unused.int2")= unused2,
         Rcpp::Named("unused.int3")= unused3,
         Rcpp::Named("non.meas.pts.flg")= non_meas_pts_flg,
         Rcpp::Named("abs.z.flg")= absZAxisFlag,
         Rcpp::Named("reserv")= reserved,
         Rcpp::Named("num.bits.pt")= num_bits_pt,
         Rcpp::Named("min.z")= minPt,
         Rcpp::Named("max.z")= maxPt,
         Rcpp::Named("num.pts.line")= num_pts_profile,
         Rcpp::Named("num.lines")= num_profiles,
         Rcpp::Named("num.pts.total")= totalNumPts
        );
    List info2;
    info2 = List::create(
         Rcpp::Named("x.inc")= x_inc,
         Rcpp::Named("y.inc")= y_inc,
         Rcpp::Named("z.inc")= z_inc,
         Rcpp::Named("x.label")= xname,
         Rcpp::Named("y.label")= yname,
         Rcpp::Named("z.label")= zname,
         Rcpp::Named("x.unit")= xstepunits,
         Rcpp::Named("y.unit")= ystepunits,
         Rcpp::Named("z.unit")= zstepunits,
         Rcpp::Named("x.len.unit")= xunits,
         Rcpp::Named("y.len.unit")= yunits,
         Rcpp::Named("z.len.unit")= zunits,
         Rcpp::Named("x.ratio")= xratio,
         Rcpp::Named("y.ratio")= yratio,
         Rcpp::Named("z.ratio")= zratio,
         Rcpp::Named("replica")= replica,
         Rcpp::Named("inverted")= inverted,
         Rcpp::Named("leveled")= leveled
        );
    List info3;
    info3 = List::create(
         Rcpp::Named("unused.sing1")= unused4,
         Rcpp::Named("unused.sing2")= unused5,
         Rcpp::Named("unused.sing3")= unused6,
         Rcpp::Named("meas.time.sec")= secs,
         Rcpp::Named("meas.time.min")= mins,
         Rcpp::Named("meas.time.hr")= hrs,
         Rcpp::Named("meas.time.day")= day,
         Rcpp::Named("meas.time.month")= month,
         Rcpp::Named("meas.time.yr")= year,
         Rcpp::Named("meas.time.wkday")= wkday,
         Rcpp::Named("meas.dur")= duration,
         Rcpp::Named("unused.int4")= unused7,
         Rcpp::Named("unused.int5")= unused8,
         Rcpp::Named("unused.sing4")= unused9,
         Rcpp::Named("unused.int6")= unused10,
         Rcpp::Named("leng.comment.zone")= leng_comment_zone,
         Rcpp::Named("leng.private.zone")= leng_private_zone,
         Rcpp::Named("free.zone")= freeZone
        );
    List info4;
    info4 = List::create(
         Rcpp::Named("x.off")= xoff,
         Rcpp::Named("y.off")= yoff,
         Rcpp::Named("z.off")= zoff,
         Rcpp::Named("T.spacing")= Tspacing,
         Rcpp::Named("T.off")= Toff,
         Rcpp::Named("T.name")= TName,
         Rcpp::Named("T.step.unit")= TUnit,
         //Rcpp::Named("comment.zone")=,
         //Rcpp::Named("private.zone")=,
         Rcpp::Named("surface") = real_surface
        );


    List infoTotal;
    infoTotal = List::create(
         Rcpp::Named("info1")= info1,
         Rcpp::Named("info2")= info2,
         Rcpp::Named("info3")= info3,
         Rcpp::Named("info4")= info4
        );

    return infoTotal;

  } catch(int e) {

    cerr<<"File open error!";

  }

}

//----------------------------------------------------
//Index to subscripts function for a rectangular matrix
//----------------------------------------------------
int * index2sub(int idx, int num_row)
{

  int r = ((idx-1) % num_row) + 1;
  int c = floor((idx-1) / num_row) + 1;
  int * subscript = new int[2];
  subscript[0] = r;
  subscript[1] = c;
  //int * scrp = &subscript;

  //return(scrp);
  return(subscript);
}
