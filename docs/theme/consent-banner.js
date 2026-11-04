silktideCookieBannerManager.updateCookieBannerConfig({
  background: {
    showBackground: false
  },
  cookieIcon: {
    position: "bottomLeft"
  },
  cookieTypes: [
    {
      id: "analytics",
      name: "Analytics",
      description: "<p>We use Google Analytics to measure how many people visit this page. This data is not used for advertising, only to give us insights into how many researchers are interested in rastair. This might help us with grant applications to justify that we're serving the community.</p>",
      required: false,
      onAccept: function() {
        gtag('consent', 'update', {
          analytics_storage: 'granted',
        });
        dataLayer.push({
          'event': 'consent_accepted_analytics',
        });
      },
      onReject: function() {
        gtag('consent', 'update', {
          analytics_storage: 'denied',
        });
      }
    }
  ],
  text: {
    banner: {
      description: "<p>We use Google analytics to measure engagement with rastair, not to serve advertisements. This helps us demonstrate that we serve the scientific community. If you would prefer not to be counted, you can safely reject cookies here.</p>",
      acceptAllButtonText: "Accept",
      acceptAllButtonAccessibleLabel: "Accept cookies",
      rejectNonEssentialButtonText: "Reject",
      rejectNonEssentialButtonAccessibleLabel: "Reject tracking",
      preferencesButtonText: "Preferences",
      preferencesButtonAccessibleLabel: "Toggle preferences"
    },
    preferences: {
      title: "Customize your cookie preferences",
      description: "<p>We respect your right to privacy. You can choose not to allow some types of cookies. Your cookie preferences will apply across our website.</p>",
      creditLinkText: "Get this banner for free",
      creditLinkAccessibleLabel: "Get this banner for free"
    }
  }
});
