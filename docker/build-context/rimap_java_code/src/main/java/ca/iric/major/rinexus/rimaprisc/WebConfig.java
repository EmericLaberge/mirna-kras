/*
 * Copyright (c) 2025 François Major, Major Lab (Université de Montréal)
 * Licensed under the MIT License. See LICENSE file in the project root for details.
 */

package ca.iric.major.rinexus.rimaprisc;

import org.springframework.context.annotation.Bean;
import org.springframework.context.annotation.Configuration;
import org.springframework.web.servlet.config.annotation.CorsRegistry;
import org.springframework.web.servlet.config.annotation.EnableWebMvc;
import org.springframework.web.servlet.config.annotation.WebMvcConfigurer;

@Configuration                                                                                                                                                          
public class WebConfig {
    @Bean
    public WebMvcConfigurer corsConfigurer() {
        return new WebMvcConfigurer() {                                                                                                                                                  
            @Override
            public void addCorsMappings( CorsRegistry registry ) {
                registry.addMapping( "/api/**" )  // Allow CORS for this endpoint
		    .allowedOrigins( "https://rimap-risc.major.iric.ca" )  // Allow only your frontend
		    .allowedMethods( "GET", "POST", "PUT", "DELETE", "OPTIONS" ) // Allow these HTTP methods
		    .allowedHeaders( "*" ) // Allow all headers
		    .allowCredentials( true ); // Allow credentials if needed
            }
        };
    }
}
