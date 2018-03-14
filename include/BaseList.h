#ifndef BASELIST_H
#define BASELIST_H

#include "TString.h"
#include <vector>
#include "AntarcticaGeometry.h" 

namespace BaseList{



  /** allows us to treat both static bases and paths the same */ 
  class abstract_base 
  {

    public: 

      virtual const char * getName() const = 0; 
      virtual const char * getSource() const = 0; 
      virtual AntarcticCoord  getPosition(unsigned time) const = 0;
      virtual bool isValid(unsigned time) const { (void) time; return true; }

      virtual ~abstract_base() { ; } 
      virtual void Draw(const char *opt = "m")  const= 0; 
  }; 


  class base : public abstract_base {
  public:
    base(const TString& theName, const TString& source, double lat, double lon, double alt=0)
      : name(theName), dataSource(source), position(AntarcticCoord::WGS84, lat, lon, alt) {;}
    base(const TString& theName, double lat, double lon, double alt=0)
      : name(theName), dataSource(""), position(AntarcticCoord::WGS84, lat, lon, alt) {;}

    virtual ~base() { ; } 

    TString name;
    TString dataSource;
    AntarcticCoord position;

    virtual AntarcticCoord getPosition(unsigned t) const {
      (void) t;
      return position.as(AntarcticCoord::WGS84);
    }
    virtual const char * getName() const { return name.Data(); } 
    virtual const char * getSource() const { return dataSource.Data(); } 
    virtual void Draw(const char * opt = "p") const; 
  };


  /** A path is a flight or traverse */ 
  class path : public abstract_base 
  {
    public:
      path(const TString & name, TString & source, 
                 int npoints, const double  * lat, const double * lon,
                 const double * alt,  const unsigned  * time) ; 
      virtual ~path() {; } 


    TString name; 
    TString dataSource;
    Bool_t isFlight; /// true for flight, false for traverse

    std::vector<AntarcticCoord> ps; 
    std::vector<unsigned> ts; 

    virtual const char * getSource() const { return dataSource.Data(); } 
    virtual const char * getName() const { return name.Data(); }
    virtual AntarcticCoord  getPosition(unsigned time) const; 
    virtual bool isValid(unsigned time) const { return time >= ts[0] && time < ts[ts.size()-1] ; } 
    virtual void Draw(const char * opt = "lp") const; 

    /** 
     * Boolian match function for std::find, checks the call signs match
     * @return true if the names match
     */
    bool operator()(const path& other){
      return name == other.name;
    }
  };


  /** Return the ith base . This function knows about the ANITA version */ 
  const base& getBase(UInt_t i);

  /** Return the ith path . This function knows about the ANITA version */ 
  const path& getPath(UInt_t i);

  /** Return the ith base or path. This function knows about the ANITA version */ 
  const abstract_base & getAbstractBase(UInt_t i); //both 


  void makeBaseList();      //refills base lists if empty. This does both paths and bases.
  void makeEmptyBaseList(); //makes base/path lists empty. Not sure why you would ever do this, but this is called somewhere... 
  
  /** Returns the number of stationary bases. ANITA-version aware */
  size_t getNumBases();

  /** Returns the number of flights/ traverses. ANITA-version aware */
  size_t getNumPaths();

  /** Returns the number of flights/ traverses + bases. ANITA-version aware */
  size_t getNumAbstractBases(); 


  /** Ever wanted to locate a base inedex? now you can! 
   * This just does a strcasestr... no regexes or anything like that sadly. 
   *
   * Returns the first match. 
   *
   * if a pointer to the a vector of indices is passed, then will fill tha twith all matching bases. 
   *
   **/ 
  int findBases(const char * query, std::vector<int> * all_matches = 0, bool include_paths = false); 

};


#endif // BASELIST_H
