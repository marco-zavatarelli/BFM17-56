#ifdef	DEBUG

#define	TRACE(format,args)		write(*,format) args 
#define	IFTRACE(condition,format,args)	if (condition) write(*,format) args 
#define	TESTNANVECTOR(args,t,s,g)		call testnan_vector(args,t,s,g)
#define	TESTNAN(scalar,gr,t,s,g)		call testnan(scalar,gr,t,s,g)
#define CHECKFLUX(A,B,C,D)              call unicflux(A,B,C,D)

#else

#define	TRACE(args)	
#define	IFTRACE(args)
#define	TESTNANVECTOR(args,t,s,g)
#define	TESTNAN(args,gr,t,s,g)
#define CHECKFLUX(A,B,C,D)

#endif


#ifdef BFM_GOTM

#define BFM_ERROR_MSG gotm_error_msg
#define BFM_ERROR gotm_error

#else

#define BFM_ERROR_MSG bfm_error_msg
#define BFM_ERROR bfm_error

#endif






