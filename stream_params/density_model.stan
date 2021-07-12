// Author: Sergey Koposov Email: skoposov _AT_ ed.ac.uk
// If you use this model, please acknowledge it by citing Erkal+2017
// https://ui.adsabs.harvard.edu/abs/2017MNRAS.470...60E/abstract
// and Koposov+2019 https://ui.adsabs.harvard.edu/abs/2019MNRAS.485.4726K/abstract papers
// History
// 12 Jul 2021 It was discovered that the spline coefficient code
// had an off-by-one bug that manifested itself when the knot
// separations were not equal. (that does not affect the Li+2021 paper
// results as we used constant separation knots)

functions{

	// get the vector of spacings between nodes
	vector geths(int n_nodes, vector nodes)
	{
		int n = n_nodes -1;
		vector[n] hs;
		for (i in 1:n)
		{
			hs[i] = nodes[i+1]-nodes[i];
		}
		return hs;
	}
	// obtain the vector of spline coefficients given the location
	// of the nodes and values there
	// We are using natural spline definition
	vector getcoeffs(int n_nodes, vector nodes, vector vals)
	{
		int n=n_nodes-1;
		vector[n] hi;
		vector[n] bi;
		vector[n-1] vi;
		vector[n-1] ui;
		vector[n_nodes] ret;
		vector[n-1] zs;
		matrix[n-1,n-1] M = rep_matrix(0, n-1, n-1);

		n = n_nodes-1;

		for (i in 1:n)
		{
			hi[i] = nodes[i+1]-nodes[i];
			bi[i] =  1/hi[i]*(vals[i+1]-vals[i]);
		}
		for (i in 2:n)
		{
			vi[i-1] = 2*(hi[i-1]+hi[i]);
			ui[i-1] = 6*(bi[i] - bi[i-1]);
		}
		for (i in 1:n-1)
		{
			M[i,i] = vi[i];
		}
		for (i in 1:n-2)
		{
			M[i+1,i] = hi[i+1];
			M[i,i+1] = hi[i+1];
		}
		//print (M)
		zs = M \ ui ; //mdivide_left_spd(M, ui);
		ret[1]=0;
		ret[n_nodes] =0;
		ret[2:n_nodes-1]=zs;

		return ret;

	}

	// Evaluate the spline, given nodes, values at the nodes
	// spline coefficients, locations of evaluation points
	// and integer bin ids of each point
	vector spline_eval(int n_nodes, vector nodes,
	       vector vals, vector zs,
	       int n_dat, vector x, int[] i)
	{

		vector[n_nodes-1] h;
		vector[n_dat] ret;
		int i1[n_dat];
		for (ii in 1:n_dat)
		{
			i1[ii] = i[ii] + 1;
		}
		h = geths(n_nodes, nodes);

		ret = (
		    zs[i1] ./ 6 ./ h[i] .* square(x-nodes[i]) .*(x-nodes[i])+
		    zs[i]  ./ 6 ./ h[i] .* square(nodes[i1]-x) .* (nodes[i1]-x)+
		    (vals[i1] ./ h[i] - h[i] .* zs[i1] ./ 6) .* (x-nodes[i])+
		    (vals[i] ./ h[i] - h[i] .* zs[i] ./ 6) .* (nodes[i1]-x)
		    );
		return ret;
	}

	// find in which node interval we should place each point of the vector
	int[] findpos(int n_nodes, vector nodes, int n_dat, vector x)
	{
		int ret[n_dat];
		for (i in 1:n_dat)
		{
			for (j in 1:n_nodes-1)
			{
				if ((x[i]>=nodes[j]) && (x[i]<nodes[j+1]))
				{
					ret[i] = j;
				}
			}
		}
		return ret;
	}
}

data {
     // number of pixels in the density map
     int n_pix;
     // number counts of objects in pixels
     int hh[n_pix];
     // the x locations of pixels
     vector[n_pix] x;
     // the y locations of pixels
     vector[n_pix] y;

     // number of nodes in the intencity spline
     int n_int_nodes;
     // location of nodes or the intensity spline
     vector[n_int_nodes] int_nodes;

     // number of nodes in the spatial density spline
     int n_fi2_nodes;
     // location of nodes for the spatial density spline
     vector[n_fi2_nodes] fi2_nodes;

     // number of nodes in the stream width spline
     int n_width_nodes;
     // location of nodes for the width spline
     vector[n_width_nodes] width_nodes;

     // number of nodes for the background density spline
     int n_bg_nodes;
     // location of nodes the background density spline
     vector[n_bg_nodes] bg_nodes;

     // the same thing as above for the background slope
     int n_bgsl_nodes;
     vector[n_bgsl_nodes] bgsl_nodes;

     // the same thing as above for the quadratic background slope
     int n_bgsl2_nodes;
     vector[n_bgsl2_nodes] bgsl2_nodes;

		 // the number of pixels we are actually fitting
     int nfit;
		 // the ids of pixels we are fitting
     int fitset[nfit];

		 // the number of pixels for which are just predicting
     int npred;
		 // the ids of pixels we are predicting
     int predset[npred];

		 // the vertical bin size
     real binsizey;
		 // the width0 of the prior
     real width_prior_cen;
		 // the width of the log width prior
     real lwidth_prior_width;
}
transformed data
{
	// These are the calculation of in which node interal each pixel
	// should go to. Since the location of nodes is different for the
	// density/track/width etc the locations ids will be different for
	// each of those parameters
	int node_ids_int[n_pix] = findpos(n_int_nodes, int_nodes, n_pix, x);
	int node_ids_fi2[n_pix] = findpos(n_fi2_nodes, fi2_nodes, n_pix, x);
	int node_ids_width[n_pix] = findpos(n_width_nodes, width_nodes, n_pix, x);
	int node_ids_bg[n_pix] = findpos(n_bg_nodes, bg_nodes, n_pix, x);
	int node_ids_bgsl[n_pix] = findpos(n_bgsl_nodes, bgsl_nodes, n_pix, x);
	int node_ids_bgsl2[n_pix] = findpos(n_bgsl2_nodes, bgsl2_nodes, n_pix, x);
	int hhfit[nfit] = hh[fitset];
	int hhpred[npred] = hh[predset];
}
parameters
{
	// The parameters of the model are the values at the
	// spline nodes for the log(intensity) of the stream,
	// log(width), log(bg), bg slope, bg quadratic slope,
	// and spatial track
	vector<lower=-100>[n_int_nodes] log_ints;
	vector[n_width_nodes] log_widths;
	vector[n_bg_nodes] log_bgs;
	vector[n_bgsl_nodes] bgsls;
	vector[n_bgsl2_nodes] bgsls2;
	vector[n_fi2_nodes] fi2s;
}

transformed parameters
{
	// The evaluation of the spline coefficients
 	// given the values at the spline nodes
	vector[n_int_nodes] coeffs_int= getcoeffs(n_int_nodes, int_nodes, log_ints);
	vector[n_fi2_nodes] coeffs_fi2= getcoeffs(n_fi2_nodes, fi2_nodes, fi2s);
	vector[n_width_nodes] coeffs_width= getcoeffs(n_width_nodes, width_nodes, log_widths);
	vector[n_bg_nodes] coeffs_bg = getcoeffs(n_bg_nodes, bg_nodes, log_bgs);
	vector[n_bgsl_nodes] coeffs_bgsl = getcoeffs(n_bgsl_nodes, bgsl_nodes, bgsls);
	vector[n_bgsl2_nodes] coeffs_bgsl2 = getcoeffs(n_bgsl2_nodes, bgsl2_nodes, bgsls2);
	vector[n_pix] logint_val;
	vector[n_pix] fi2_val;
	vector[n_pix] logwidth_val;
	vector[n_pix] logbg_val;
	vector[n_pix] bgsl_val;
	vector[n_pix] bgsl2_val;

	vector[n_pix] width_val;
	vector[n_pix] xmod;
	vector[n_pix] logbg_pix;
	vector[n_pix] logint_pix;

	vector[n_width_nodes] widths;

	// Actual evaluation of the splines at each pixel
	logint_val = spline_eval(n_int_nodes, int_nodes,
	       log_ints, coeffs_int,
	       n_pix, x, node_ids_int);
	fi2_val = spline_eval(n_fi2_nodes, fi2_nodes,
	       fi2s, coeffs_fi2,
	       n_pix, x, node_ids_fi2);
	logwidth_val = spline_eval(n_width_nodes, width_nodes,
	       log_widths, coeffs_width,
	       n_pix, x, node_ids_width);
	logbg_val = spline_eval(n_bg_nodes, bg_nodes,
	       log_bgs, coeffs_bg,
	       n_pix, x, node_ids_bg);
	bgsl_val = spline_eval(n_bgsl_nodes, bgsl_nodes,
	       bgsls, coeffs_bgsl,
	       n_pix, x, node_ids_bgsl);
	bgsl2_val = spline_eval(n_bgsl2_nodes, bgsl2_nodes,
	       bgsls2, coeffs_bgsl2,
	       n_pix, x, node_ids_bgsl2);

	// log densities of the background/stream at each pixel
	logbg_pix = logbg_val + bgsl_val/10 .* y + bgsl2_val/100 .* y .* y;
	//logint_pix = logint_val - 0.5 * square(y-fi2_val) ./ exp(2 * logwidth_val);
//	logint_pix = logint_val  + log(normal_cdf(y+0.05, fi2_val, exp(logwidth_val))-
//		     		 normal_cdf(y-0.05,fi2_val,exp(logwidth_val))+1e-100);
	width_val = exp(logwidth_val);
	for (i in 1:n_pix)
	{

		logint_pix[i] = logint_val[i] + log ( width_val[i] *sqrt(2*pi()))  + log(
							normal_cdf(y[i] + binsizey/2., fi2_val[i], width_val[i])-
		     		 	normal_cdf(y[i] - binsizey/2., fi2_val[i], width_val[i])
							+1e-100);
	}
        //print (log(normal_cdf(0+0.05, 1, 0.05)-
        //                         normal_cdf(0-0.05,1,0.05)+1e-100));
	// - 0.5 * square(y-fi2_val) ./ exp(2 * logwidth_val);
	for (i in 1:n_pix)
	{
	//	print(i," ",logint_pix[i]," ",y[i]," ",fi2_val[i]);
		xmod[i] = log_sum_exp(logbg_pix[i],logint_pix[i]);
	}
	widths = exp(log_widths);

}
model
{
	// Priors
	target += normal_lpdf(fi2s|0, 1);
	target += normal_lpdf(log_widths|log(width_prior_cen),lwidth_prior_width);

	//Likelihood
	hhfit ~ poisson_log(xmod[fitset]);
}

generated quantities
{
	real log_lik;
	real log_lik_pred;
	// Save the actuall log-likelihood of the function evaluation
	log_lik = poisson_log_lpmf(hh|xmod);
	log_lik_pred = poisson_log_lpmf(hh[predset]|xmod[predset]);
}
